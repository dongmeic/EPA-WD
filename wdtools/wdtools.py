# TODO: replace paths with pathlib.Path
from collections import Counter
import difflib
import io
import json
import os
from os import walk
import pickle
import re
import string
import time
import warnings
import webbrowser

import fiona
import geopandas as gpd
from geopy.geocoders import Nominatim
import googlemaps
import numpy as np
import openpyxl
import pandas as pd
from PyPDF2 import PdfReader, PdfWriter
import requests
#from shapely.geometry import shape
from shapely.validation import make_valid
from shapely.errors import ShapelyDeprecationWarning
from win32com.client import Dispatch

from const import (
    ALL_TXID, COL_DICT, COUNTY_DICT, INPATH, OR_COUNTIES, OUTPATH, TAXLOT_PATH,
    TRANSFORMER, TRSQQ, TRSQQ_DICT, TSQ_DST, TTDF, VAR_LIST, WD_PATH)
from deliverables import SASplitter
from taxlots import MapIndexReader, TaxlotReader
from utils import (
    check_duplicates, create_ORMap_name, get_lot_numbers, reindex_data,
    remove_duplicates, split_WD_to_taxmaps)


warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)


# Clean up for now; will deal with placement later...
google_key = json.load(open('config/keys.json'))['google_maps']['APIKEY']
#outfolder = 'test'
outfolder = 'output\\matched'
trsqq_correction_dict = dict(
    zip(
        list(range(0, 6)),
        ['township number', 'township direction', 'range number',
         'range direction', 'section number', 'QQ']))
nm2add = [0, 1420, 2143, 2878, 3932, 4370]
issuepath = f'{INPATH}\\GIS\\ArcGIS Pro Project\\DataReview\\issueIDs.gdb'


def read_trsqq():
    'Read trsqq list, dictionary, and dataframe'    
    with open(os.path.join(INPATH, 'trsqq_list.pkl'), 'rb') as f:
        trsqq = pickle.load(f)
    with open(os.path.join(INPATH, 'trsqq_dict.pkl'), 'rb') as f:
        trsqq_dict = pickle.load(f) 
    df = pd.read_csv(os.path.join(INPATH, 'trsqq_df.csv'))
    return trsqq, trsqq_dict, df

# More to deal with later...
cnts = gpd.read_file(f'{INPATH}\\GIS\\Oregon_Counties.shp')
pdf_outpath = r'L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping\output\pdf'

with open(os.path.join(INPATH, 'ParticipCnt.pkl'), 'rb') as f:
    county_list = pickle.load(f)


# Taxlot Review -----------------------------------------------------------
def read_taxlots(year):
    return TaxlotReader(county_list).read(year)

readTaxlots = read_taxlots  # for backwards compatibility


def read_map_index(year):
    return MapIndexReader().read(year)

ReadMapIndex = read_map_index


# Deliverable -------------------------------------------------------------
def format_gdf_provided(gdf, wdID):
    'Format gdf provided to be in the same format as mapped records'
    gdf['wdID'] = wdID
    gdf = gdf.dissolve('wdID')
    gdf['wdID'] = gdf.index
    gdf.reset_index(drop=True, inplace=True)
    gdf['code'] = 3
    gdf = gdf.to_crs(epsg=2992)
    if 'Shape_Length' not in gdf.columns:
        gdf['Shape_Length'] = gdf.length
    if 'Shape_Area' not in gdf.columns:
        gdf['Shape_Area'] = gdf.area
    gdf = gdf[['wdID', 'code', 'Shape_Length', 'Shape_Area', 'geometry']]
    return gdf


def get_corrected_wd_df(num, setID=False):
    'Combine all the corrected DSL tables into one dataframe'
    if setID:
        wd_df = pd.read_csv(
            fr'{WD_PATH}\Corrected_by_Set\Set{num}_2017-20220601.csv')
    else:
        frames = []
        for i in range(num):
            wd_dt = pd.read_csv(
                fr'{WD_PATH}\Corrected_by_Set\Set{i + 1}_2017-20220601.csv')
            frames.append(wd_dt)
        wd_df = pd.concat(frames, ignore_index=True)
    return wd_df


def export_wd_gdf_by_record(gdf, out_name):
    'Export geodataframe with the original wd data'
    gdf = gdf[VAR_LIST]
    gdf['received_date'] = gdf['received_date'].dt.date
    gdf['response_date'] = gdf['response_date'].dt.date
    gdf['lat'], gdf['lon'] = TRANSFORMER.transform(
        gdf.centroid.x, gdf.centroid.y)
    gdf = gdf.rename(columns=COL_DICT)
    try:
        gdf.to_file(fr'{OUTPATH}\test\{out_name}.shp')
    except RuntimeError:
        gdf['geometry'] = gdf.geometry.buffer(0)
    gdf.to_file(fr'{OUTPATH}\test\{out_name}.shp')
    return gdf


def split_sa_by_rid(
        wd_df, sa_gdf_all, all_map_idx, all_taxlot, em_wids, do_export=False,
        out_name='example_data', do_review=False):
    return SASplitter(
        VAR_LIST, COL_DICT, OUTPATH, wd_df, sa_gdf_all, all_map_idx,
        all_taxlot, em_wids, do_export=False, out_name='example_data',
        do_review=False
    ).split_by_rd()


split_SA_by_rid_in_df = split_sa_by_rid


def read_all_map_idx(do_export=False):
    'Combine mapIndex from all year'
    frames = []
    for year in range(YEAR_START, YEAR_END):
        map_idx_dt = read_map_idx(year)
        map_idx_dt['year'] = str(year)
        frames.append(map_idx_dt[['year', 'ORMapNum', 'geometry']])
    df = pd.concat(frames, ignore_index=True)
    gdf = gpd.GeoDataFrame(df, crs='EPSG:2992', geometry='geometry')
    if do_export:
        with open(os.path.join(INPATH, 'ORMapIndex.pkl'), 'wb') as f:
            pickle.dump(list(gdf.ORMapNum.unique()), f) 
    return gdf
    

def read_map_idx(year):
    'Read mapIndex from one year'
    map_idx = gpd.read_file(
        fr'{INPATH}\GIS\ORMAP_data\ORMAP_Taxlot_Years\Taxlots{year}.gdb', 
        layer='MapIndex')
    return map_idx


def replace_geometry(gdf):
    '''Replace the geometry of SA ploygons with manual review
    Args:
    - gdf (geopandas.GeoDataFrame): geodataframe containing the SA polygons to
      replace
    (rid is record ID)
    '''
    path = f'{INPATH}\GIS\ArcGIS Pro Project\DataReview\DataReview.gdb'
    new_p_layers = [layer for layer in fiona.listlayers(path) if 'rid' in layer]
    frames = []
    for new_p_layer in new_p_layers:
        new_p_layer = gpd.read_file(path, layer=new_p_layer)
        rid = int(new_p_layer.replace('rid', ''))
        new_p_layer['record_ID'] = rid
        frames.append(new_p_layer[['record_ID', 'geometry']])
    sa_df = pd.concat(frames, ignore_index=True)
    sa_gdf = gpd.GeoDataFrame(sa_df, geometry='geometry')
    unique_rids = sa_df.record_ID.unique()
    df = gdf.drop(columns=['geometry'])
    gdf1 = df[df.record_ID.isin(unique_rids)].merge(sa_gdf,on='record_ID')
    gdf2 = pd.concat(
        [gdf[~gdf.record_ID.isin(unique_rids)], gdf1], ignore_index=True)
    return gdf2


def get_all_wd(n, is_raw=False):
    'Combine all the original DSL tables into one dataframe'    
    frames = []
    for i in range(n):
        wd_dt = combine_wd_tables(
            setID=f'Set00{i + 1}', n_to_add=nm2add[i], is_raw=is_raw)
        wd_dt['SetID'] = i + 1
        frames.append(wd_dt)
    wd_df = pd.concat(frames, ignore_index=True)
    return wd_df


## HERE: move to utils
def combine_wd_tables(setID, n_to_add, is_raw=True):
    '''Combine all the wd tables in the set to review unique records, record_ID is
    used for combined tables.
    Use this function when reindex is not neccessary
    Args:
    -n_to_add (int) the number of previous records (records from the previous sets;
        same for all functions with this variable)
    '''
    if is_raw:
        frames = []
        files = list_files(os.path.join(WD_PATH, setID))
        # in case there are unidentified files
        files = [f for f in files if '~$' not in f]
        for f in files:
            data_file = os.path.join(WD_PATH, setID, f)
            xl = pd.ExcelFile(data_file)
            wd_dt = pd.read_excel(data_file, sheet_name=xl.sheet_names[1])
            frames.append(wd_dt)
        wd_df = pd.concat(frames, ignore_index=True)
    else:
        wd_df = pd.read_csv(fr'{WD_PATH}\Corrected_by_Set\{setID}.csv')
    # this creates unique IDs for all the records in the same set
    wd_df['record_ID'] = range(1, wd_df.shape[0] + 1) 
    wd_df['record_ID'] = wd_df.copy().record_ID + n_to_add
    selected_ids = wd_df.parcel_id.astype(str) != 'nan'
    wd_df.loc[selected_ids, 'notes'] = (
        wd_df[selected_ids]['parcel_id'].apply(lambda x: make_notes(x)))
    # get year from the receive date
    if raw:
        wd_df.loc[:, 'recyear'] = wd_df.received_date.apply(lambda x: x.year)
    else:
        wd_df.loc[:, 'recyear'] = wd_df.received_date.apply(
            lambda x: int(x.split('-')[0]))
    # get year from the wd ID 
    wd_df['IDyear'] = wd_df.wetdet_delin_number.apply(lambda x: x[2:6]) 
    selected_ids = wd_df.parcel_id.astype(str) != 'nan'
    wd_df.loc[selected_ids, 'missinglot'] = wd_df[selected_ids].parcel_id.apply(
        lambda x: without_lots(x))
    return wd_df


def make_notes(text):
    """
    add notes to the wd records
    """
    r = re.search('ROW|RR', text, re.IGNORECASE)
    p = re.search('partial|part|p|portion', text, re.IGNORECASE)
    m = re.search('Many|multiple|SEVERAL|various', text, re.IGNORECASE)
    w = re.search('Water', text, re.IGNORECASE)
    l = re.search('RAIL', text, re.IGNORECASE)
    res = []
    for srch, msg in zip([r, p, m, w, l], ['ROW','Partial','Many','Water','Rail']):
        if srch:
            #res += msg
            res.append(msg)
    res = ', '.join(res)
    return res


def without_lots(text):
    """
    check whether the parcel ID is without lots
    """
    if any(c.isdigit() for c in text):
        nms = re.findall(r'\d+', text)
        if all([any([x in text for x in [nm+'st', nm+'nd', nm+'rd', nm+'th',
                                                  nm+' st', nm+' th', nm+' nd', nm+' rd']]) for nm in nms]):
            res = 'Y'
        else:
            res = 'N'
    else:
        res = 'Y'
    return res



def get_added_SA():
    path = f'{INPATH}\\GIS\\ArcGIS Pro Project\\DataReview\\added.gdb'
    lyrs = [lyr for lyr in fiona.listlayers(path)]
    frames = []
    for lyr in lyrs:
        sa = gpd.read_file(path, layer=lyr)
        sa = sa.to_crs(epsg=2992)
        frames.append(sa)
    sadf = pd.concat(frames, ignore_index=True)
    return sadf

def get_all_SA(num):
    """
    combine all the SA polygons into one dataframe
    allsa is to check all SA polygons after reviewing issue IDs
    """  
    frames = []
    for i in range(num):
        file1 = OUTPATH + f'\\final\\mapped_wd_Set00{i+1}.shp'
        file2 = OUTPATH + f'\\final\\Set00{i+1}_mapped_wd.shp'
        if os.path.exists(file2):
            sa_dt = gpd.read_file(file2)
        else:
            sa_dt = gpd.read_file(file1)
        sa_dt['SetID'] = i+1
        frames.append(sa_dt)
    sa_df = pd.concat(frames, ignore_index=True)
    sa_gdf = gpd.GeoDataFrame(sa_df, geometry='geometry')
    added = get_added_SA()
    sa_gdf = pd.concat([sa_gdf[['wdID', 'code', 'geometry']], 
                   added[['wdID', 'code', 'geometry']]], 
                   ignore_index=True)
    return sa_gdf

def join_WD_with_SA_by_taxmap(df, gdf, mapindex, outnm='wd_mapped_data', export=True):
    """
    df is the corrected WD dataframe from get_all_wd
    gdf is the SA polygons from get_all_SA
    mapindex is the combined taxmaps from read_all_mapIdx
    """
    frames = []
    wdlist = gdf.wdID.unique()
    df['ORMapNum'] = df[['county', 'trsqq']].apply(lambda row: create_ORMap_name(county=row.county, trsqq=row.trsqq), axis = 1)
    for wid in wdlist:
        #print(wid)
        df_s = df[df.wetdet_delin_number==wid]
        gdf_s = gdf[gdf.wdID==wid]
        ORmn = df_s.ORMapNum.values
        yr = wid[2:6]
        mapIdx = mapindex[(mapindex.ORMapNum.isin(ORmn)) & (mapindex.year==yr)]
        exp_gdf = split_WD_to_taxmaps(df=df_s, gdf=gdf_s, wdID=wid, mapindex=mapIdx)
        exp_gdf['lat'], exp_gdf['lon'] = TRANSFORMER.transform(exp_gdf.representative_point().x, exp_gdf.representative_point().y)
        ndf = df_s.drop_duplicates(subset='ORMapNum')
        ndf.drop(columns=['parcel_id','site_id','record_ID'], inplace=True)
        g = df_s.groupby('ORMapNum')
        pi_df = g.apply(lambda x: '; '.join(x.parcel_id.astype(str).unique())).to_frame(name='parcel_id').reset_index(drop=True)
        si_df = g.apply(lambda x: '; '.join(x.site_id.astype(str).unique())).to_frame(name='site_id').reset_index(drop=True)
        ri_df = g.apply(lambda x: '; '.join(x.record_ID.astype(str).unique())).to_frame(name='record_ID').reset_index()
        sdf = pd.concat([ri_df, pi_df, si_df], axis=1)
        fdf = ndf.merge(sdf, on='ORMapNum')
        exp_gdf = fdf.merge(exp_gdf[['code','ORMapNum','lat','lon','geometry']], on='ORMapNum')
        frames.append(exp_gdf)
    sa_df = pd.concat(frames, ignore_index=True)
    sa_gdf = gpd.GeoDataFrame(sa_df, geometry='geometry')
    if export:
        sa_gdf=sa_gdf.rename(columns=COL_DICT)
        try:
            sa_gdf.to_file(os.path.join(INPATH, "output", "final", f"{outnm}.shp"), index=False)
        except RuntimeError:
            sa_gdf['geometry'] = sa_gdf.geometry.buffer(0)
            sa_gdf.to_file(os.path.join(INPATH, "output", "final", f"{outnm}.shp"), index=False)
    return sa_gdf
 
################################################ Report #########################################################

def flatten(l):
    """
    convert lists in list to a list
    """
    return [item for sublist in l for item in sublist]

def count_lst_ele(dct, lst, ctnm):
    """
    count list elements
    ctnm - count number
    """
    qaqc_cnt = [*map(dct.get, lst)]
    freq = Counter(qaqc_cnt)
    df = pd.DataFrame(sorted(freq.items()))
    df.columns = ['county', ctnm]
    return df

def read_list(setid):
    """
    read a list 
    """
    with open(os.path.join(INPATH, f"{setid}_mapped.pkl"), "rb") as f:
        lst = pickle.load(f)
    return lst

# clean all formats
def removeFormatting(ws):
    """
    remove Excel format
    """
    # ws is not the worksheet name, but the worksheet object
    for row in ws.iter_rows():
        for cell in row:
            cell.style = 'Normal'
    return ws

def reformat(file):
    """
    reformat Excel
    """    
    wb = openpyxl.load_workbook(file)
    ws = wb.active
    ws = removeFormatting(ws)
    wb.save(file)
    print("Removed format")
    excel = Dispatch('Excel.Application')
    wb = excel.Workbooks.Open(file)
    excel.Worksheets(1).Activate()
    excel.ActiveSheet.Columns.AutoFit()
    wb.Close(True)
    print("Autofitted columns...")
    
################################################ Tier 3 & 4 #####################################################

def writelist(lst, lstnm, setID):
    """
    write a list to a pickle file
    """
    with open(os.path.join(INPATH, f'{setID}_{lstnm}.pkl'), 'wb') as f:
            pickle.dump(lst, f)
    
def review_loop_r1(setID=None, wdid_list=None, df=None, partial=False, idx=False, wd_id=None, wddf=None):
    """
    loop through the unmatched records and check the original records
    """
    toadd=[]
    if not partial:
        df = pd.read_csv(os.path.join(
            INPATH, f'\\output\\to_review\\unmatched_df_{setID}_r1_N.csv'))
    if df is not None:
        wdid_list = list(df.wetdet_delin_number.unique())
    else:
        if wdid_list is None:
            print("need to provide a list of wdID..")
            return None
    n = len(wdid_list)
    if wd_id is None:
        i = -1
    else:
        i =  wdid_list.index(wd_id)
    for wdID in wdid_list[i+1:]:
        j = wdid_list.index(wdID)
        print(f'{round(((j/n)*100),1)}% digitized, {n-j} records remained, expected to be done in about {int(((n-j)*0.8)+0.5)} hours...')
        print(wdID)
        if idx:
            if df is not None:
                print(f'index = {df[df.wetdet_delin_number==wdID].index[0]+1}')
                print(check_unmatched_r1(wdID = wdID, df = df))
            else:
                print(f'index = {wdid_list.index(wdID)+1}')
                print(check_unmatched_r1(wdID = wdID, df = wddf)) 
        user_input = input("Press 'p' to pause or any key to stop...")
        if user_input in ['p', 'P']:
            while True:
                user_input = input("Press 'a' to add the wd record or 'c' to continue...")
                if user_input in ['a', 'A']:
                    toadd.append(wdID)
                    break
                if user_input in ['c', 'C']:
                    break
        else:
            break
        time.sleep(1)
    if toadd != []:
        return toadd, wdID

def check_unmatched_r1(wdID, df):
    """
    check unmatched records for a given WDID
    df is from the splited unmatched records from r1 process
    """
    url = df.loc[df.wetdet_delin_number == wdID, 'DecisionLink'].values[0]
    selcols = ['county', 'trsqq', 'parcel_id', 'latitude', 'longitude', 'record_ID', 'notes', 'missinglot', "status_name", "is_batch_file"]
    if str(url) == 'nan':
        print('Decision link is not available')
    else:
        webbrowser.open(url)
    return df.loc[df.wetdet_delin_number == wdID, selcols]

def review_loop(setID):
    """
    loop through the unmatched records and check the review notes
    df is from r2 notes
    """
    df = pd.read_csv(
        os.path.join(
            INPATH,
            '\\output\\to_review\\',
            f'review_unmatched_{setID}_r2_N_0.csv'))
    df = df.reset_index()
    for i in range(df.shape[0]):
        wdID = df.loc[i, 'wetdet_delin_number']
        print(wdID)
        print(check_review_notes_r2n(wdID = wdID, df = df))
        user_input = input("Press 'p' to pause or any other key to continue...")
        if user_input in ['p', 'P']:
            while True:
                user_input = input("Press 'c' to continue...")
                if user_input in ['c', 'C']:
                    break
        time.sleep(1)

def check_review_notes_r2n(wdID, df):
    """
    check review notes for a given WDID
    df is from r2 notes
    """
    url = df.loc[df.wetdet_delin_number == wdID, 'DecisionLink'].values[0]
    if str(url) == 'nan':
        print('Decision link is not available')
    else:
        webbrowser.open(url)
    return df.loc[df.wetdet_delin_number == wdID, ['correct_type', 'correction', 'cor_trsqq']]

def check_completeness(setID='003', a=3):
    """
    check completeness of the mapping
    a: numbers to add, for the records that are edited in the matched records directly
    """
    revpath = f'{INPATH}\GIS\ArcGIS Pro Project\DataReview\Set{setID}.gdb'
    partial = gpd.read_file(revpath, layer=f'{setID}_partial')
    mapped1 = list(partial.wdID.unique())
    mapped0 = [lyr for lyr in fiona.listlayers(revpath) if (lyr not in [f'Set{setID}_wo_lot', f'Set{setID}_partial']) and ('L' not in lyr)]
    mapped2 = list(map(lambda x: x.replace('_', '-'), mapped0))
    mapped = mapped1 + mapped2
    matched = gpd.read_file(
        f'{INPATH}\\output\matched\matched_records_Set{setID}_edited.shp')
    pct = (len(mapped)+a)/len(matched[~matched.notes.isnull()].wdID.unique())
    print(f'{round(pct*100, 1)}% completed...')
    return sorted(mapped), len(mapped)

def extract_page_from_locPath(filePath, pageNm, wdID, k=0):
    """
    Extract a page from a pdf file from a local path
    """
    pdf_file = PdfReader(filePath)
    pageObj = pdf_file.getPage(pageNm-1)
    pdf_writer = PdfWriter()
    pdf_writer.addPage(pageObj)
    output = f'{pdf_outpath}\\{wdID}_{pageNm+k}.pdf'
    with open(output, 'wb') as output_pdf:
        pdf_writer.write(output_pdf) 

def extract_page_from_docLink(url, pageNm, wdID):
    """
    Extract a page from a pdf file from a url
    """
    if str(url) == 'nan':
        print('Decision link is not available')
    else:
        response = requests.get(url=url, timeout=120)
        on_fly_mem_obj = io.BytesIO(response.content)
        pdf_file = PdfReader(on_fly_mem_obj)
        pageObj = pdf_file.getPage(pageNm-1)
        pdf_writer = PdfWriter()
        pdf_writer.addPage(pageObj)
        output = f'{pdf_outpath}\\{wdID}_{pageNm}.pdf'
        with open(output, 'wb') as output_pdf:
            pdf_writer.write(output_pdf) 

def review_mapped(setID):
    """
    Review the mapped taxlots
    """
    revpath = f'{INPATH}\GIS\ArcGIS Pro Project\DataReview\{setID}.gdb'
    mapped0 = [lyr for lyr in fiona.listlayers(revpath) if (lyr not in [f'{setID}_wo_lot', f'{setID}_partial']) and ('L' not in lyr)]
    for wID in mapped0:
        gdf = gpd.read_file(revpath, layer=wID)
        if 'wdID' in gdf.columns:
            if len(gdf.wdID.unique()) > 1:
                print(wID)

def rename_wdID(x):
    """
    Replace underscore with dash
    """
    if '_' in x:
        x = x.replace('_', '-')
    return x

def revise_single_partial_file(setID, wID, from_set=True):
    """
    Revise the geometry of a single partial taxlot
    """
    if from_set:
        revpath = f'{INPATH}\GIS\ArcGIS Pro Project\DataReview\{setID}.gdb'
    else:
        revpath = issuepath   
    gdf = gpd.read_file(revpath, layer=wID)
    gdf = gdf.to_crs(epsg=2992)
    selcols = ['Shape_Length', 'Shape_Area']
    if 'wdID' in gdf.columns:
        df = gdf.copy()[selcols +  ['wdID']]
        df.loc[:, 'wdID'] = df.wdID.apply(lambda x: rename_wdID(x))
    else:
        df = gdf.copy()[selcols]
        df.loc[:,'wdID'] = wID.replace('_', '-')
    df.loc[:,'geometry'] = gdf.loc[:,'geometry']
    return df

def merge_single_partial_file(setID, wIDlist, from_set=True):
    """
    Merge the revised partial taxlots into a single GeoDataFrame
    """
    df = pd.DataFrame()
    if from_set:
        for wID in wIDlist:
            df=pd.concat([df, revise_single_partial_file(setID, wID)], ignore_index=True)
    else:
        for wID in wIDlist:
            df=pd.concat([df, revise_single_partial_file(setID=False, wID=wID, from_set=False)], ignore_index=True)
            
    gdf = gpd.GeoDataFrame(df, crs="EPSG:2992", geometry='geometry')
    return gdf

def combine_matched_digitized(setID, editedIDs, nm_to_add, export=True):
    """
    Combine the edited matched records, digitized partial taxlots, taxlots without lot IDs, and the list of issue IDs
    """
    revpath = f'{INPATH}\GIS\ArcGIS Pro Project\DataReview\{setID}.gdb'
    # get separated feature files
    mapped0 = [lyr for lyr in fiona.listlayers(revpath) if (lyr not in [f'{setID}_wo_lot', f'{setID}_partial']) and ('L' not in lyr)]
    matched = gpd.read_file(
        f'{INPATH}\\output\matched\matched_records_{setID}_edited.shp')
    # get edited feature in the original matches
    edited_gdf = matched[matched.wdID.isin(editedIDs)]
    edited_gdf = edited_gdf[['wdID', 'geometry']].dissolve('wdID')
    edited_gdf.loc[:, 'wdID'] = edited_gdf.index
    edited_gdf['code'] = 1
    # get digitized partially-matched files
    partial = gpd.read_file(revpath, layer=f'{setID}_partial')
    partial = partial.to_crs(epsg=2992)
    gdf = merge_single_partial_file(setID=setID, wIDlist=mapped0)
    dat = gdf.append(partial, ignore_index=True)
    if dat.shape[0] > len(dat.wdID.unique()):
        dat = dat[['wdID', 'geometry']].dissolve('wdID')
        dat.loc[:, 'wdID'] = dat.index
    # merge the ones without lot IDs
    wo_lot = gpd.read_file(revpath, layer=f'{setID}_wo_lot')
    wo_lot = wo_lot.to_crs(epsg=2992)
    # merge digitized or edited features
    data1 = wo_lot.append(dat[['wdID', 'geometry']], ignore_index=True)
    data1['code'] = 2
    data2 = edited_gdf.append(data1[['wdID', 'code', 'geometry']], ignore_index=True)
    # exclude the ones that were digitized or edited in the original matches
    excluded = [wdID for wdID in data2.wdID.unique() if wdID in matched.wdID.unique()]
    issues = pd.read_csv(os.path.join(INPATH, "output", "to_review", f"{setID}_Mapping_Issues.csv"))
    # exclude the ones that have issues
    issueIDs = list(issues.wetdet_delin_number.unique())
    file = OUTPATH+f"\\matched\\{setID}_reviewed.txt"
    if os.path.exists(file):
        with open(file) as f:
            reviewedIDs = f.readlines()
        issueIDs = [iID for iID in issueIDs if iID not in reviewedIDs]
    file = OUTPATH+f"\\matched\\{setID}_not_mapped.txt"
    if os.path.exists(file):
        with open(file) as f:
            withdrawnIDs = f.readlines()
            data2 = data2[~data2.wdID.isin(withdrawnIDs)]
        issueIDs = issueIDs + withdrawnIDs
    matched_gdf = matched[~matched.wdID.isin(excluded+issueIDs)]
    matched_gdf['code'] = 0
    # merge matched and digitized/edited
    data3 = matched_gdf[['wdID', 'code', 'geometry']].append(data2[['wdID', 'code', 'geometry']], ignore_index=True)
    shp = data3.copy()
    shp.geometry = shp.apply(lambda row: make_valid(row.geometry) if not row.geometry.is_valid else row.geometry, axis=1)
    final_gdf = shp.dissolve('wdID')
    final_gdf.loc[:, 'wdID'] = final_gdf.index
    wd = combine_wd_tables(setID=setID, n_to_add=nm_to_add)
    unmatchedIDs = [wdID for wdID in wd.wetdet_delin_number.unique() if wdID not in final_gdf.wdID.unique()]
    toCheck = [ID for ID in unmatchedIDs if ID not in issueIDs]
    digitized_nIDs = len(editedIDs) + len(dat.wdID.unique()) + len(wo_lot.wdID.unique())
    file = OUTPATH+f"\\matched\\{setID}_edited_1.txt"
    if os.path.exists(file):
        with open(file) as f:
            editedIDs1 = f.readlines()
        final_gdf.loc[final_gdf.wdID.isin(editedIDs1), 'code'] = 1
    if export:
        final_gdf.to_file(os.path.join(
            INPATH, "output", "final", f"mapped_wd_{setID}.shp"), index=False)
    return final_gdf, toCheck, matched_gdf, digitized_nIDs, unmatchedIDs, issueIDs

def run_Tier3_4_final(setID, nm_to_add):
    """
    combine matched and digitized records and review
    gdf: the final shapefile that combined both automatic matches and digitized records
    """
    start = time.time()
    with open(OUTPATH+f'\\matched\\{setID}_edited.txt') as f:
        edited = f.readlines()
    file = OUTPATH+f"\\matched\\{setID}_edited_1.txt"
    if os.path.exists(file):
        with open(file) as f:
            editedIDs1 = f.readlines()
    edited = edited + editedIDs1
    gdf, toCheck, matched_gdf, digitized_nIDs, unmatchedIDs, issueIDs = combine_matched_digitized(setID=setID, 
                                                                                     editedIDs=edited[0].split(", "), 
                                                                                     nm_to_add=nm_to_add)
    end = time.time()
    print(f'it took {round((end - start)/60, 0)} minutes to complete')
    return gdf, toCheck, matched_gdf, digitized_nIDs, unmatchedIDs, issueIDs

################################################ Tier 2 #####################################################
def get_point_from_lonlat(lon, lat, transprj=True, export=True):
    """
    Get a point from a lon/lat pair
    """
    df = pd.DataFrame([[lon, lat]], columns=['Longitude', 'Latitude'])
    gdf = gpd.GeoDataFrame(df, crs="EPSG:4326", geometry=gpd.points_from_xy(df.Longitude, df.Latitude))
    if transprj:
        gdf = gdf.to_crs(epsg=2992)
    if export:
        gdf.to_file(f'{INPATH}\\test\point.shp')
    return gdf

# point in polygon - WD point in taxtlot
# require taxtlot
def extract_taxlot_info(wd_pt, taxlot, year):
    """
    Extract taxlot information from a point in polygon analysis
    """
    pip_mask = taxlot.contains(wd_pt.loc[0, 'geometry'])
    if any(pip_mask):
        pip_data = taxlot.loc[pip_mask].copy()
    else:
        pip_data = gpd.sjoin_nearest(wd_pt, taxlot, distance_col='dist')
    pip_data.loc[:, 'YDiff'] = abs(pip_data.year.astype(int) - int(year))
    ID = pip_data.loc[pip_data.YDiff == np.min(pip_data.YDiff.values), 'ORTaxlot'].values[0]
    if not any(pip_mask):
        dist = pip_data.loc[pip_data.YDiff == np.min(pip_data.YDiff.values), 'dist'].values[0]
        ID = ID + f', about {int(dist+0.5)} ft away'
    return ID

def search_for_county_name(reverse_geocode_result):
    my_list = reverse_geocode_result[0]['address_components']
    search_key = 'long_name'
    search_value = 'County'
    for my_dict in my_list:
        if search_key in my_dict and search_value in my_dict[search_key]:
            return my_dict[search_key].split(' '+search_value)[0]
            break
    else:
        print('County was not found on Maps!')
        return None 

def extract_county_name(point, polygon):
    """
    Extract county name from a point in polygon analysis
    """
    pip_mask = polygon.contains(point.loc[0, 'geometry'])
    if any(pip_mask):
        pip_data = polygon.loc[pip_mask].copy()
    else:
        pip_data = gpd.sjoin_nearest(point, polygon, distance_col='dist')
    cnty = pip_data.NAME.values[0]
    if not any(pip_mask):
        dist = pip_data.dist.values[0]
        print(f'Found {cnty}, about {int(dist+0.5)} degree away')
    return cnty    
    
def get_county_code_from_lonlat(lon, lat, search="OR"):
    """
    Get county code from longitude and latitude.
    The search method includes "GM" (GoogleMaps), "OSM" (OpenStreetMaps), and "OR" (Shapefile for Counties in Oregon, default)
    """
    if search == "GM":
        gmaps = googlemaps.Client(key=google_key)
        reverse_geocode_result = gmaps.reverse_geocode((lat, lon))
        cnty = search_for_county_name(reverse_geocode_result)
    elif search == "OSM":   
        geolocator = Nominatim(user_agent="geoapiExercises")
        location = geolocator.reverse(str(lat)+","+str(lon))
        address = location.raw['address']
        county = address.get('county', '')
        cnty = county.split(' ')[0]
    else:
        pnt = get_point_from_lonlat(lon, lat, transprj=False, export=False)
        cnty = extract_county_name(pnt, cnts)
    res = str(COUNTY_DICT[cnty]).zfill(2)
    return res

def separate_numbers_letters(x):
    """
    Separate numbers and letters in a string.
    """
    numbers = re.findall('\d+', x)
    letters = re.findall("[a-zA-Z]+", x)
    return numbers, letters

def find_different_indices(list1, list2):
    """
    Find the indices of different items in two lists.
    """
    different_indices = []
    for i, (a, b) in enumerate(zip(list1, list2)):
        if a != b:
            different_indices.append(i)
    return different_indices

def pad_string(string, length=10):
    """
    Pad a string with 0s to a specified length.
    """
    if len(string) < length:
        padded_string = string.ljust(length, '0')
    else:
        padded_string = string
    return padded_string

def combine_lists(list1, list2):
    combined_list = list(zip(list1, list2))
    return combined_list

def remove_tuple_format(input_list):
    """
    Remove tuple format from a list.
    """
    output_list = []
    for item in input_list:
        if isinstance(item, tuple):
            output_list.extend(item)
        else:
            output_list.append(item)
    return output_list

def get_list_elements_by_index(input_list, index_list):
    """
    Get elements from a list by index.
    """
    output_list = [input_list[i] for i in index_list]
    return output_list

def split_trsqq(trsqq_to_check):
    """
    Split a TRSQQ into a list of numbers and letters.
    """
    numbers1, letters1 = separate_numbers_letters(trsqq_to_check[:-2])
    letters1.append(trsqq_to_check[-2:])
    trsqq_to_check_lst = remove_tuple_format(combine_lists(numbers1, letters1))
    return trsqq_to_check_lst

def compare_trsqq(trsqq_to_check, trsqq_to_compare):
    """
    Compare two TRSQQs and return the indices of the different elements, the correct elements, and the errors.
    """
    trsqq_to_check_lst = split_trsqq(trsqq_to_check)
    trsqq_to_compare_lst = split_trsqq(trsqq_to_compare)
    diff_idx = find_different_indices(trsqq_to_check_lst, trsqq_to_compare_lst)
    if len(diff_idx) != 0:
        if len(diff_idx) == 1:
            correct_trsqq_elements = trsqq_to_compare_lst[diff_idx[0]]
            errors = trsqq_to_check_lst[diff_idx[0]]
        else:
            correct_trsqq_elements = get_list_elements_by_index(trsqq_to_compare_lst, diff_idx)
            errors = get_list_elements_by_index(trsqq_to_check_lst, diff_idx)
    return diff_idx, correct_trsqq_elements, errors

def join_list_elements(my_list):
    """
    function to join list elements
    """
    delimiter1 = ', '
    delimiter2 = ' and '
    position1 = 1
    position2 = -1
    res = delimiter1.join(my_list[:position1]) + delimiter1 + delimiter1.join(my_list[position1:position2]) + delimiter2 + my_list[position2]
    return res

def report_trsqq_correction(trsqq_to_check, trsqq_to_compare, to_correct=False):
    """
    function to review and correct trsqq
    """
    diff_idx, correct_trsqq_elements, errors = compare_trsqq(trsqq_to_check, trsqq_to_compare)
    if len(diff_idx) == 1:
        keylist = trsqq_correction_dict[diff_idx[0]]
        res = f'from {errors} to {correct_trsqq_elements}'
        if to_correct:
            return trsqq_to_compare
        else:
            return keylist, res
    else:
        keylist = [trsqq_correction_dict.get(key) for key in diff_idx]
        if len(keylist) == 2:
            joined_keys = ' and '.join(keylist)
            joined_errors = ' and '.join(errors)
            joined_corrections = ' and '.join(correct_trsqq_elements)
            res = f'from {joined_errors} to {joined_corrections}'
        else:
            joined_keys = join_list_elements(keylist)
            res = f'from {join_list_elements(errors)} to {join_list_elements(correct_trsqq_elements)}'
        if to_correct:
            return trsqq_to_compare
        else:
            return joined_keys, res

def has_letter(string):
    """
    Checks if a string contains any letter.
    Returns True if the string contains at least one letter, False otherwise.
    """
    return any(char.isalpha() for char in string)

def get_lot_number_from_taxlot(x):
    if 'ROAD' in x:
        res = 'ROADS'
    elif 'WATER' in x:
        res = 'WATER'
    elif 'RAIL' in x:
        res = 'RAILS'
    else:
        lot = x.split('--')[1]
        if has_letter(lot):
            res = lot
        else:
            res = str(int(lot))
    return res

def correct_trsqq(trsqq_to_check, lon, lat, taxlot, year):
    """
    Corrects trsqq based on the taxlot information derived from the coordinate.
    """
    #print(trsqq_to_check)
    wd_pt  = get_point_from_lonlat(lon = lon, lat = lat)
    tID = extract_taxlot_info(wd_pt = wd_pt, taxlot = taxlot, year = year)
    trsqq_to_check_c = pad_string(trsqq_to_check)
    trsqq_to_compare = taxlot2trsqq(tID)
    trsqq_to_compare_c = pad_string(trsqq_to_compare)
    res = report_trsqq_correction(trsqq_to_check_c, trsqq_to_compare_c, to_correct=True)    
    return res, tID

# this function only works when the coordinates are accurate and one-on-one match among WD ID, trsqq, parcel IDs, and coordinate
# limitation - one WD record (possibly with multiple records with different trasqq and parcel IDs) has only one coordindate
def review_wd_record_w_coord(wd_id, county_to_check, trsqq_to_check, parcel_IDs_to_check, lon, lat, taxlot, year):
    """
    Reviews a WD record with coordinate information.
    input includes all the key fields to be checked for taxlot information
    """
    print(f'reviewing {wd_id}')
    wd_pt  = get_point_from_lonlat(lon = lon, lat = lat)
    tID = extract_taxlot_info(wd_pt = wd_pt, taxlot = taxlot, year = year)
    if "away" not in tID:
        trsqq_to_check_c = pad_string(trsqq_to_check)
        lots_to_check = get_lot_numbers(parcel_IDs_to_check)
        trsqq_to_compare = taxlot2trsqq(tID)
        trsqq_to_compare_c = pad_string(trsqq_to_compare)
        lots_to_compare = TTDF.loc[TTDF.trsqq==trsqq_to_compare, 'ORTaxlot'].values
        lots_to_compare = list(map(get_lot_number_from_taxlot, lots_to_compare))
        if trsqq_to_compare_c == trsqq_to_check_c:
            print("trsqq matched, checking county code...")
            cnty_code = int(get_county_code_from_lonlat(lon, lat))
            county_to_compare = [key for key, value in COUNTY_DICT.items() if value == cnty_code][0]
            # need to check the typos in the county name first
            if county_to_check == county_to_compare:
                print("county code is corrected, need to check lot numbers...")
                if any([x not in lots_to_compare for x in lots_to_check]):
                    lots_to_correct = [x for x in lots_to_check if x not in lots_to_compare]
                    cor_type, cor_notes = "lot number", f'lot number {lots_to_correct} might be incorrect, the matched taxlot is {tID} for {trsqq_to_compare}'
                    print("lot numbers might be wrong...")
                else:
                    notes = 'the record seems to be correct, need to review why it was not matched in match_wd_data_with_taxlot'
                    print(notes)
                    cor_type, cor_notes = None, notes
            else:
                cor_type, cor_notes = "county", f'from {county_to_check} to {county_to_compare}'
                print("corrected county...")
        else:
            # check the lot numbers
            lots_matched = [x for x in lots_to_check if x in lots_to_compare]
            if len(lots_matched) > 0:
                if len(lots_matched) == len(lots_to_check):
                    print("all lots are matched...")
                    cor_type, cor_notes = report_trsqq_correction(trsqq_to_check_c, trsqq_to_compare_c)
                    print("corrected trsqq...")
                else:
                    notes = f"some lots are not matched, need to review trsqq, the close-match is {trsqq_to_compare}"
                    print(notes)
                    print(f"lots to check: {lots_to_check}, and lots to compare: {lots_to_compare}")
                    cor_type, cor_notes = "to review", notes
            else:
                notes = f"there is not any matched lot, need to review trsqq, the close-match is {trsqq_to_compare}"
                print(notes)
                print(f"lots to check: {lots_to_check}, and lots to compare: {lots_to_compare}")
                cor_type, cor_notes = "to review", notes
    else:
        notes = f'coordinate might be incorrect, nearby taxlot is {tID}'
        print(notes)
        cor_type, cor_notes = "coordinate", notes
    return cor_type, cor_notes

def split_unmatched_df(df, ml, setID, export=True):
    """
    split unmatched df into two parts: one part is the records with multiple matches, the other part is the records with only one match
    df: the output from report_unmatched
    ml (missing lot): whether the unmatched records are missing parcel id in digit
    setID: the set ID of the unmatched records
    """
    df = df[df.missinglot == ml]
    IDcol = 'wetdet_delin_number'
    value_counts = df[IDcol].value_counts()
    r1_df = df[df[IDcol].isin(value_counts[value_counts > 1].index)]
    r2_df = df[df[IDcol].isin(value_counts[value_counts == 1].index)]
    r2_df = r2_df[((r2_df.latitude.astype(str) != 'nan')|(r2_df.longitude.astype(str) != 'nan'))]
    if export:
        r1_df.to_csv(
            os.path.join(
                INPATH,
                '\\output\\to_review\\',
                f'unmatched_df_{setID}_r1_{ml}.csv'),
            index=False)
        r2_df.to_csv(os.path.join(INPATH, '\\output\\to_review\\', f'unmatched_df_{setID}_r2_{ml}.csv'), index=False)
    return r1_df, r2_df

def generate_taxlot_output(df, taxlot, setID, ml, export=False):
    """
    generate the taxlot output for the unmatched records
    df: the second output from split_unmatched_df
    taxlot: the taxlot shapefile
    setID: the setID of the unmatched records
    ml (missing lot): whether the unmatched records are missing parcel id in digit
    """
    df['pairs'] = list(zip(df['IDyear'].astype(str), df['ORTaxlot']))
    taxlots_to_review = taxlot[taxlot[['year', 'ORTaxlot']].apply(tuple, axis=1).isin(df.pairs.values)]
    taxlots_to_review_2 = taxlots_to_review.merge(df, on='ORTaxlot')
    taxlots_to_review_2.drop(columns=['pairs'], inplace=True)
    if export:
        taxlots_to_review_2.to_file(os.path.join(
            INPATH,
            '\\output\\to_review\\',
            f'review_unmatched_{setID}_r2_{ml}.shp'))
    return taxlots_to_review_2

def taxlot_from_coord(x):
    """
    get trsqq and txid from the correction notes
    """
    txid = x.split('the matched taxlot is ')[1].split(' for ')[0]
    trsqq = x.split(' for ')[1]
    return trsqq, txid

def trsqq_from_nearby_taxlot(x):
    """
    get trsqq and txid from the correction notes
    """
    txid = x.split(', about')[0].replace('coordinate might be incorrect, nearby taxlot is ', '')
    trsqq = taxlot2trsqq(txid)
    return trsqq, txid

def review_unmatched_df_r2(df, taxlot, setID, ml, export=True):
    """
    df: the second output from split_unmatched_df
    taxlot: the taxlot shapefile
    setID: the ID of the set
    ml (missing lot): whether the unmatched records are missing parcel id in digit
    """
    # exclude the unusual county records, e.g., Yamhill and Washington
    df = df[df.county.isin(OR_COUNTIES)]
    outdf = df.copy()[['wetdet_delin_number', 'trsqq', 'parcel_id', 'county', 'latitude', 'longitude', 'DecisionLink', 'record_ID', 'IDyear']]
    outdf.loc[:,'correct_type'], outdf.loc[:,'correction'] = zip(*outdf.apply(lambda row: review_wd_record_w_coord(wd_id = row.wetdet_delin_number, 
                                                                        county_to_check = row.county, 
                                                                        trsqq_to_check = row.trsqq, 
                                                                        parcel_IDs_to_check = row.parcel_id, 
                                                                        lon = row.longitude, 
                                                                        lat = row.latitude, 
                                                                        taxlot = taxlot, 
                                                                        year = row.IDyear), axis = 1))
    sel = ~outdf.correct_type.isin(['county', 'lot number', 'coordinate', None])
    outdf.loc[sel, 'cor_trsqq'], outdf.loc[sel, 'ORTaxlot'] = zip(*outdf.loc[sel,:].apply(lambda row: correct_trsqq(trsqq_to_check = row.trsqq, 
                                                                               lon = row.longitude,
                                                                               lat = row.latitude,
                                                                               taxlot = taxlot,
                                                                               year = row.IDyear), axis = 1))
    if 'lot number' in outdf.correct_type.unique():
        sel1 = outdf.correct_type=='lot number'
        outdf.loc[sel1, 'cor_trsqq'], outdf.loc[sel1, 'ORTaxlot'] = zip(*outdf.loc[sel1, 'correction'].apply(lambda x: taxlot_from_coord(x)))
    
    if 'coordinate' in outdf.correct_type.unique():
        sel2 = outdf.correct_type=='coordinate'
        outdf.loc[sel2, 'cor_trsqq'], outdf.loc[sel2, 'ORTaxlot'] = zip(*outdf.loc[sel2, 'correction'].apply(lambda x: trsqq_from_nearby_taxlot(x)))

    if export:
        outdf.to_csv(
            os.path.join(
                INPATH,
                '\\output\\to_review\\',
                f'review_unmatched_{setID}_r2_{ml}_0.csv'),
            index=False)
    return outdf

trsqq_cor_dict = {'township number': 0, 
                  'township direction':2, 
                  'range number':3, 
                  'range direction':5,
                  'section number':6, 
                  'QQ':8}

def replace_str_index(text,index=0,replacement=''):
    """
    replace a character in a string at a given index
    """
    return text[:index] + replacement + text[index+len(replacement):]

def adjust_cor(x):
    """
    adjust the correction from single digit to double digit
    """
    if str(x).isdigit():
        return str(x).zfill(2)
    else:
        return x

def adjust_cor_from_to(df):
    """
    adjust the correction data frame where the from and to columns are single digit to double digit
    df: reindexed df from reindex_data; df from correct_unmatched or update_unmatched_df_r2 or combine_corrected_unmatched function
    """
    df.loc[:, 'from'] = df['from'].apply(lambda x: adjust_cor(x))
    df.loc[:, 'to'] = df['to'].apply(lambda x: adjust_cor(x))
    return df

# need to run review_unmatched_df_r2 and do some manual review work first to get the notes
# the document to review is review_unmatched_Set00?_r2_N_0.csv
def correct_unmatched(df, setID, s, ml, export=True):
    """
    df: the output from split_unmatched_df
    setID: the set ID
    s: the split number - r1 or r2, r for review
    ml (missing lot): whether the unmatched records are missing parcel id in digit
    export: whether to export the corrected unmatched records
    """
    notes = pd.read_csv(
        os.path.join(
            INPATH,
            '\\output\\to_review\\',
            f'unmatched_df_{setID}_{s}_{ml}_notes.csv'))
    notes = adjust_cor_from_to(notes)
    rID = 'record_ID' in notes.columns
    df = df.copy()[df.wetdet_delin_number.isin(notes.wetdet_delin_number.unique())]
    for wdID in df.wetdet_delin_number.unique():
        
        if rID:
            for r_id in notes[notes.wetdet_delin_number == wdID]['record_ID'].unique():
                sel = (df.wetdet_delin_number == wdID) & (df.record_ID == r_id)
                sel1 = (notes.wetdet_delin_number == wdID) & (notes.record_ID == r_id)
                fields = notes.loc[sel1, 'field'].unique()
                for field in fields:
                    sel2 = sel1 & (notes.field==field)
                    if field == 'trsqq':
                        df.loc[sel, field] = df.loc[sel, field].apply(lambda x: pad_string(x))
                        cor_types = notes.loc[sel2, 'cor_type'].values
                        for cor_type in cor_types:
                            sel3 = sel2 & (notes.cor_type==cor_type)
                            repval = notes.loc[sel3, 'to'].values[0]
                            valrep = notes.loc[sel3, 'from'].values[0]
                            if ('.5' in str(valrep)) or (cor_type == 'trsqq'):
                                df.loc[sel, field] = df.loc[sel, field].apply(lambda x: x.replace(valrep, repval))
                            else:
                                ind = trsqq_cor_dict[cor_type]
                                df.loc[sel, field] = df.loc[sel, field].apply(lambda x: replace_str_index(x, index=ind, replacement=repval))          
                    else:
                        valrep = notes.loc[sel2, 'from'].values[0]
                        repval = notes.loc[sel2, 'to'].values[0]
                        df.loc[sel, field] = df.loc[sel, field].apply(lambda x: x.replace(valrep, repval))
        else:    
            sel = df.wetdet_delin_number == wdID
            fields = notes.loc[notes.wetdet_delin_number == wdID, 'field'].unique()
            for field in fields:
                sel2 = (notes.wetdet_delin_number==wdID) & (notes.field==field)
                if field == 'trsqq':
                    df.loc[sel, field] = df.loc[sel, field].apply(lambda x: pad_string(x))
                    cor_types = notes.loc[sel2, 'cor_type'].values
                    for cor_type in cor_types:
                        sel3 = (notes.wetdet_delin_number==wdID) & (notes.field==field) & (notes.cor_type==cor_type)
                        ind = trsqq_cor_dict[cor_type]
                        df.loc[sel, field] = df.loc[sel, field].apply(lambda x: replace_str_index(x,
                                                                                            index=ind,
                                                                                            replacement=notes.loc[sel3, 'to'].values[0]))
                else:
                    to = notes.loc[sel2, 'to'].values[0]
                    # if to.isdigit():
                    #     to = to.zfill(2)
                    df.loc[sel, field] = df.loc[sel, field].apply(lambda x: x.replace(notes.loc[sel2, 'from'].values[0], to))
    if export:
        df.to_csv(
            os.path.join(
                INPATH,
                '\\output\\to_review\\',
                f'review_unmatched_{setID}_{s}_{ml}_1.csv'),
            index=False)
    return df

# need to run review_unmatched_df_r2 and do some manual work frist to get the notes
def update_unmatched_df_r2(df, setID, ml, export=True):
    """
    combine corrected and uncorrected in the unmatched r2 records
    df: the second output from split_unmatched_df
    setID: the setID
    ml (missing lot): whether the unmatched records are missing parcel id in digit
    export: whether to export the updated unmatched df
    """
    # rev_df is the output from review_unmatched_df_r2, nt_df is from manual input
    rev_df = pd.read_csv(
        os.path.join(
            INPATH,
            '\\output\\to_review\\',
            f'review_unmatched_{setID}_r2_{ml}_0.csv'))
    nt_df = pd.read_csv(
        os.path.join(
            INPATH,
            '\\output\\to_review\\',
            f'unmatched_df_{setID}_r2_{ml}_notes.csv'))
    df = df.copy()[~df.wetdet_delin_number.isin(nt_df.wetdet_delin_number.unique())]
    rev_df = rev_df.copy()[~rev_df.wetdet_delin_number.isin(nt_df.wetdet_delin_number.unique())]
    rev_df = rev_df[['wetdet_delin_number', 'cor_trsqq']]
    ndf = df.merge(rev_df, on='wetdet_delin_number')
    selectedID = ndf.cor_trsqq.astype(str) != 'nan'
    ndf.loc[selectedID, 'trsqq'] = ndf.loc[selectedID, 'cor_trsqq'].apply(lambda x: x.rstrip('0'))
    ndf.drop(columns='cor_trsqq', inplace=True)
    rev_df = pd.read_csv(
        os.path.join(
            INPATH,
            '\\output\\to_review\\',
            f'review_unmatched_{setID}_r2_{ml}_1.csv'))
    rdf = ndf.append(rev_df, ignore_index = True)
    if export:
        rdf.to_csv(
            os.path.join(
                INPATH,
                '\\output\\to_review\\',
                f'review_unmatched_{setID}_r2_{ml}_2.csv'),
            index=False)
    return rdf

def combine_corrected_unmatched(setID, ml, skip=True, export=True):
    """
    combine the output from correct_unmatched and update_unmatched_df_r2
    setID: the setID of the unmatched records
    ml (missing lot): whether the unmatched records are missing parcel id in digit
    """
    rev_df1 = pd.read_csv(
        os.path.join(
            INPATH,
            '\\output\\to_review\\',
            f'review_unmatched_{setID}_r1_{ml}_1.csv'))
    # if skip update_unmatched_df_r2
    if skip:
        rev_df2 = pd.read_csv(
            os.path.join(
                INPATH,
                '\\output\\to_review\\',
                f'review_unmatched_{setID}_r2_{ml}_1.csv'))
    else:    
        rev_df2 = pd.read_csv(
            os.path.join(
                INPATH,
                '\\output\\to_review\\',
                f'review_unmatched_{setID}_r2_{ml}_2.csv'))
    df = rev_df1.append(rev_df2, ignore_index = True)
    if export:
        df.to_csv(
            os.path.join(
                INPATH,
                '\\output\\to_review\\',
                f'review_unmatched_{setID}_{ml}.csv'), index=False)
    return df

def review_WD_record_via_Pro(gdf, wdID):
    """
    review the unmatched records via Pro
    gdf: the output from split_unmatched_df
    wdID: the wetdet_delin_number of the unmatched records
    return the gdf of the unmatched records
    """
    gdf = gdf[gdf.wdID == wdID]
    gdf.to_file(os.path.join(INPATH, 'output', 'wd_shp', wdID + '.shp'))
    return gdf

# make sure ORTaxlot is in the right format from df
def get_taxlot_to_check_r2(revdf, taxlot, setID, ml):
    """
    revdf: the output from review_unmatched_df_r2
    taxlot: the taxlot shapefile
    setID: the setID of the unmatched records
    ml (missing lot): whether the unmatched records are missing parcel id in digit
    """
    df = revdf.copy()
    df['pairs'] = list(zip(df['IDyear'].astype(str), df['ORTaxlot']))
    taxlots_to_review = taxlot[taxlot[['year', 'ORTaxlot']].apply(tuple, axis=1).isin(df.pairs.values)]
    taxlots_to_review_2 = taxlots_to_review.merge(df, on='ORTaxlot')
    taxlots_to_review_2.drop(columns=['pairs'], inplace=True)
    taxlots_to_review_2.rename(columns={'wetdet_delin_number': 'wdID', 
                      'DecisionLink':'doc_link',
                      'correct_type':'cor_type'}, inplace=True)
    taxlots_to_review_2.to_file(
        os.path.join(
            INPATH,
            'output',
            'to_review',
            f'review_unmatched_{setID}_r2_{ml}.shp'))
    return taxlots_to_review_2

def adjust_taxlot(tx, ty):
    """
    adjust the taxlot with correct trsqq and taxlot with sheet number
    """
    res = ty
    if tx in TSQ_DST:
        pattern = ty.split('--')[1][1:]
        tid_list = TTDF[TTDF.trsqq == tx].ORTaxlot.unique()
        for tid in tid_list:
            if re.search(pattern, tid):
                res = tid
                break
    return res

def adjust_taxlot_df(df):
    """
    adjust the taxlot with correct trsqq and taxlot with sheet number in dataframe
    """
    df.loc[:, 'ORTaxlot'] = df.copy()[['trsqq', 'ORTaxlot']].apply(lambda row: adjust_taxlot(row.trsqq, row.ORTaxlot), axis=1)
    return df

def run_Tier2_step1(setID, unmatched_df, all_taxlot):
    """
    split unmatched records
    """
    r1_df, r2_df = split_unmatched_df(unmatched_df, ml='N', setID=setID)
    return r1_df, r2_df

def run_Tier2_step3(r1_df, r2_df, setID, nm_to_add, wd, all_taxlot):
    """
    update the match
    """
    df = combine_corrected_unmatched(setID, ml='N')
    rev_df = reindex_data(df)
    matched = match_wd_data_with_taxlot(rev_df, setID, all_taxlot, export=True, update=True)
    unmatched_df = report_unmatched(matched, setID, nm_to_add, mute = False)
    matched_toReview = matched[matched.notes.notnull()] 
    wd_toReview = wd[wd.wetdet_delin_number.isin(matched_toReview.wdID.unique())]
    wd_toReview.to_csv(OUTPATH + f'\\to_review\\partial_matched_{setID}.csv', index=False)
    unmatched_df.to_csv(
        os.path.join(
            INPATH,
            '\\output\\to_review\\',
            f'unmatched_df_{setID}_2.csv'),
        index=False)
    return matched, unmatched_df

def report2DSL(setID):
    """
    generate the correction report for DSL
    """
    r1_notes = pd.read_csv(
        os.path.join(
            INPATH,
            '\\output\\to_review\\',
            f'unmatched_df_{setID}_r1_N_notes.csv'))
    r2_notes = pd.read_csv(
        os.path.join(
            INPATH,
            '\\output\\to_review\\',
            f'unmatched_df_{setID}_r2_N_notes.csv'))
    cordf = r1_notes.append(r2_notes, ignore_index=True)
    matched = gpd.read_file(OUTPATH + f'\\matched\\matched_records_{setID}.shp')
    corrected = matched[matched.record_ID.isin(cordf.record_ID.values)][['wdID','trsqq', 'parcel_id', 'record_ID']].drop_duplicates(ignore_index=True)
    report = cordf.merge(corrected[['trsqq', 'parcel_id', 'record_ID']], on='record_ID')
    report['trsqq'] = report.trsqq.apply(lambda x: x.rstrip('0'))
    report.to_csv(OUTPATH + f'\\corrected\\corrected_{setID}.csv', index=False)
    return report
    
################################################ Tier 1 #####################################################
# gdf below generally refers to the matched records
def list_files(path, folder=False):
    """
    This function takes a path and returns a list of files in the path
    """
    f = []
    for (dirpath, dirnames, filenames) in walk(path):
        f.extend(filenames)
        break
    if folder:
        f = [x[0] for x in os.walk(path)]
    return f

#def unique(list1):
#    """
#    This function takes a list and returns a list of unique values
#    """
#    x = np.array(list1)
#    return list(np.unique(x))


def read_wd_table(setID, file):
    """
    read wd tables, recordID is used for single tables
    """
    datafile = os.path.join(WD_PATH, setID, file)
    xl = pd.ExcelFile(datafile)
    wd_dt = pd.read_excel(datafile, sheet_name=xl.sheet_names[1])
    wd_dt.loc[:, 'county'] = wd_dt.county.apply(lambda x: x.capitalize())
    # this will show the records in the original single table and the ID will be updated when the tables are combined
    wd_dt.loc[:, 'recordID'] = range(1, wd_dt.shape[0] + 1)
    return wd_dt

def scan_trsqq(x):
    """
    scan trsqq to find out problematic records
    """
    nms = re.findall('\d+', x)
    lts = re.findall("[a-zA-Z]+", x)
    nmerr = [len(nm) > 2 for nm in nms[0:2]]
    lterr = [lt not in ["E", "W", "S", "N"] for lt in lts[0:2]]
    if((len(nms) <= 2)|(len(lts) <= 2)) & any(lterr):
        print(f"trsqq {x} has less than three numbers or letters")
    elif any(nmerr) & any(lterr):
        print(f"trsqq {x} has {sum(nmerr)} wrong number(s) and {sum(lterr)} wrong letter(s)")
    elif any(nmerr):
        print(f"trsqq {x} has {sum(nmerr)} wrong number(s)")
    elif any(lterr):
        print(f"trsqq {x} has {sum(lterr)} wrong letter(s)")
    else:
        #print(f"trsqq {x} seems correct")
        res=0
        return res
    

def clean_wd_table(setID, file):
    """
    function to clean up wd data by single file
    """
    wd_dt = read_wd_table(setID, file)
    # this will help identify problematic records with numbers
    selectedID = wd_dt.parcel_id.astype(str) != 'nan'
    wd_dt.loc[selectedID, 'notes'] = wd_dt.copy()[selectedID]['parcel_id'].apply(lambda x: make_notes(x))
    wd_dt['county'] = wd_dt['county'].apply(lambda x: x.title())
    wd_dt = wd_dt[wd_dt.county.isin(OR_COUNTIES)]
    if wd_dt.empty:
        res = wd_dt
    else:
        ndf = reindex_data(wd_dt)
        # get year from the receive date
        ndf.loc[:, 'recyear'] = ndf.copy().received_date.apply(lambda x: x.year)
        # get year from the wd ID 
        ndf.loc[:, 'IDyear'] = ndf.wetdet_delin_number.apply(lambda x: x[2:6])
        ndf.loc[selectedID, 'missinglot'] = ndf.loc[selectedID, 'parcel_id'].apply(lambda x: without_lots(x))
        ndf['response_date'] = ndf['response_date'].dt.strftime("%Y-%m-%d")
        ndf['received_date'] = ndf['received_date'].dt.strftime("%Y-%m-%d")
        ndf['county'] = ndf['county'].apply(lambda x: x.title())
        ndf = ndf[ndf.county.isin(OR_COUNTIES)]
        res = ndf
        #print(f'cleaned up wd data in {file} and it took about {end - start} seconds')
    return res


def read_taxlot(year, mute=True):
    """
    function to read taxlot data
    """
    txfilepath = os.path.join(TAXLOT_PATH, 'Taxlots' + str(year) + '.gdb')
    start = time.time()
    tx_dt = gpd.read_file(txfilepath, layer='TL_Dissolv')
    end = time.time()
    if not mute:
        print(f'got taxlot data in {year} and it took about {str(round((end - start)/60, 0))} minutes')
    return tx_dt

def merge_data_by_year(setID, file, year):
    """
    merge wd table with taxlot polygons in a single year
    use this function when checking data by year
    """
    wd_dt = clean_wd_table(setID, file)
    tx_dt = read_taxlot(year)
    merged = wd_dt[wd_dt.IDyear == str(year)].merge(tx_dt[['ORTaxlot', 'geometry']], 
                                                    on='ORTaxlot', 
                                                    how='left')
    #print('got merged data between wd and taxlot')
    return merged

def get_record_dict(setID, wd_df):
    """
    get a county and record in a dictionary to update recordID
    input wd_df is the output from combine_wd_tables (read all set data without merging)
    """
    counties = wd_df.county.unique()
    count_records = []
    for cnty in counties:
        count = len(wd_df[wd_df.county == cnty])
        count_records.append(count)
    record_df = pd.DataFrame({'county': counties, 'rcount':count_records})
    record_df['cum_count'] = record_df.copy()[['rcount']].cumsum(axis = 0, skipna = True).rcount.values
    record_dict = dict(zip(record_df.county[1:len(counties)], record_df.cum_count[0:(len(counties)-1)]))
    return counties, record_dict

def combine_taxlot(exportID=False, 
                   YEAR_START=2016, 
                   YEAR_END=2023):
    """
    combine taxlots from all years
    """
    frames = []
    for year in range(YEAR_START, YEAR_END):
        tx_dt = read_taxlot(year)
        tx_dt['year'] = str(year)
        frames.append(tx_dt[['year', 'ORTaxlot', 'geometry']])
    df = pd.concat(frames, ignore_index=True)
    gdf = gpd.GeoDataFrame(df, crs="EPSG:2992", geometry='geometry')
    if exportID:
        with open(os.path.join(INPATH, "ORTaxlot.pkl"), "wb") as f:
            pickle.dump(list(gdf.ORTaxlot.unique()), f) 
    return gdf

def update_recordID(df, wd_df, setID, nm_to_add):
    """
    update record_ID if needed
    input 1) wd_df is the output from combine_wd_tables (read all files in the same set without merging with taxlots); 
          2) df is the reindexed wd data, from combined_reindexed_data
    """
    counties, record_dict = get_record_dict(setID, wd_df)
    selected_cnty = df.county.isin(counties[1:])
    df = df.copy()
    df.loc[selected_cnty, 'record_ID'] = df.loc[selected_cnty, 'recordID'] + df[selected_cnty].county.map(record_dict) + nm_to_add
    df.loc[df.county == counties[0], 'record_ID'] = df.loc[df.county == counties[0], 'recordID'] + nm_to_add
    df.loc[:, 'record_ID'] = df.record_ID.astype('int64', copy=False)
    return df

def combined_reindexed_data(setID, nm_to_add):
    """
    combine reindexed wd data in the same set folder
    """
    frames = []
    files = list_files(os.path.join(WD_PATH, setID))
    # in case there are unidentified files
    files = [file for file in files if '~$' not in file]
    for file in files:
        wd_dt = clean_wd_table(setID, file = file)  
        frames.append(wd_dt)
    df = pd.concat(frames, ignore_index=True)
    wd_df = combine_wd_tables(setID, nm_to_add)
    ndf = update_recordID(df, wd_df, setID, nm_to_add)
    return ndf 

def match_wd_data_with_taxlot(df, setID, all_taxlot, export=False, update=False):
    """
    get the geometry from adjacent years with the same taxlot ID
    input df is reindexed from combined_reindexed_data or reindex_data
    return geodata with matched geometry
    run this to include only the matched records with the original ID
    make sure the matched_records_{setID}.shp is not the updated version from a previous run when update is true
    """
    tocheck_txid = df.ORTaxlot.unique()
    found = [txid for txid in tocheck_txid if txid in ALL_TXID]
    # unfound = [txid for txid in tocheck_txid if txid not in all_txid]
    # if len(unfound) > 0:
    #     sdf = df[df.ORTaxlot.isin(unfound)]
    #     sdf = adjust_taxlot_df(sdf)
    #     if sdf.shape[0] > 0:
    #         adjusted = [txid for txid in list(sdf.ORTaxlot.unique()) if txid not in unfound]
    #         rIDs_toadj = sdf[sdf.ORTaxlot.isin(adjusted)].record_ID.unique()
    #         df.loc[df.record_ID.isin(rIDs_toadj), 'ORTaxlot'] = sdf.copy()[sdf.record_ID.isin(rIDs_toadj)].ORTaxlot.values
    #         found = found + adjusted
    if len(found) > 0:
        tocheck_df = df[df.ORTaxlot.isin(found)]
        taxlot_tocheck = all_taxlot[all_taxlot.ORTaxlot.isin(found)]
        taxlot_tocheck = taxlot_tocheck.merge(tocheck_df, on='ORTaxlot', how='left')
        taxlot_tocheck.loc[:, 'ydiff'] = taxlot_tocheck.copy()[['year', 'IDyear']].apply(lambda row: abs(int(row.year) - int(row.IDyear)), axis=1)
        tdf = taxlot_tocheck.sort_values(by=['ORTaxlot', 'ydiff'])
        # keep the taxlot from the closest year
        tdf = tdf.drop_duplicates(subset='ORTaxlot', keep="first")
        ndf = tocheck_df.merge(tdf[['ORTaxlot', 'geometry']], on='ORTaxlot', how='left')
        ndf.rename(columns={'wetdet_delin_number': 'wdID', 
                      'address_location_desc':'loc_desc', 
                      'Coord-Source': 'coord_src',
                      'DocumentName':'doc_name',
                      'DecisionLink':'doc_link',
                      'is_batch_file':'isbatfile',
                      'status_name': 'status_nm',
                      'received_date':'receiveddt', 
                      'response_date':'responsedt',
                      'reissuance_response_date':'reissuance' 
                      }, inplace=True)
        ngdf = gpd.GeoDataFrame(ndf, crs="EPSG:2992", geometry='geometry')
        selcols = ['wdID', 'trsqq', 'parcel_id', 'notes', 'lots', 'lot', 'ORTaxlot', 'record_ID', 'geometry']
        if update:
            matched = gpd.read_file(
                os.path.join(
                    INPATH,
                    f'\\{outfolder}\\',
                    f'matched_records_{setID}.shp'))
            ngdf = matched.append(ngdf[selcols], ignore_index = True)
        if export: 
            # the update will overwrite the first output
            ngdf[~ngdf.geometry.isnull()][selcols].to_file(
                os.path.join(
                    INPATH,
                    f'\\{outfolder}\\',
                    f'matched_records_{setID}.shp'), 
                driver='ESRI Shapefile')  
        return ngdf
    else:
        print('no matched records found')
        return None

def report_unmatched(gdf, setID, nm_to_add, mute = True, export=False):
    """
    input gdf from match_wd_data_with_taxlot
    return the unmatched data after matching the taxlot data
    """
    wd_df = combine_wd_tables(setID, nm_to_add)
    matched_rID = gdf.record_ID.unique()
    unmatched_wd_df = wd_df[~wd_df.record_ID.isin(matched_rID)]
    df = unmatched_wd_df
    IDcol = 'wetdet_delin_number'
    counts = df[IDcol].value_counts()
    counts_df = pd.DataFrame({'ind': counts.index, 'order':range(len(counts))})
    cnt_dict = counts_df.set_index('ind').to_dict(orient='dict')['order']
    sorted_df = df.sort_values(by=[IDcol], key=lambda x: x.map(cnt_dict))
    nc = unmatched_wd_df[unmatched_wd_df.parcel_id.astype(str) == 'nan'].shape[0]
    nr = round((nc/wd_df.shape[0]) * 100, 2)
    r = round((unmatched_wd_df.shape[0]/wd_df.shape[0]) * 100, 2)
    if not mute:
        print(f'it is about {r}% of data in the original {wd_df.shape[0]} records unmatched')
        print(f'there are {nc} records ({nr}% of the original records) without parcel id')
    if export:
        sorted_df.to_csv(
            os.path.join(
                INPATH,
                '\\output\\to_review\\',
                f'unmatched_df_{setID}.csv'),
            index=False)  
    return sorted_df
    
def compare_data_report(gdf, setID, nm_to_add, export = False):
    """
    compare the output data with the existing output from the manual process
    input gdf from match_wd_data_with_taxlot
    return 1) missed_match_ID, missed_gdf: what is missing in the existing matched data; 
           2) addedID, added_gdf: what is added in the existing matched data;
           3) missed_ID: unmatched IDs
    """
    unmatched_wd_df = report_unmatched(gdf, setID, nm_to_add)
    # unmatched IDs from the run
    missed_ID = unmatched_wd_df.record_ID.unique()
    setgdf = gpd.read_file(
        os.path.join(
            INPATH, 'GIS', 'Join_Statewide.gdb'), layer=f'WD_{setID}_Combined')
    setgdf.loc[setgdf.Record_ID.astype(str) != 'nan', 'Record_ID'] = setgdf[setgdf.Record_ID.astype(str) != 'nan'].Record_ID.astype('int64', copy=False)
    matched_rID = gdf.record_ID.unique()
    # missed IDs in the existing data that is not nan
    missed_gdf = setgdf[setgdf.Record_ID.astype(str) != 'nan'][~setgdf.Record_ID.isin(matched_rID)]
    # matched IDs in the existing data that is not nan
    matched_gdf = setgdf[(setgdf.Record_ID.astype(str) != 'nan') & (setgdf.Record_ID.isin(matched_rID))]
    matchedID = matched_gdf.Record_ID.unique()
    missed_match_ID = [ID for ID in list(matched_rID) if ID not in list(matchedID)]
    missed_gdf = gdf[gdf.record_ID.isin(missed_match_ID)]
    addedID = [ID for ID in list(setgdf.Record_ID.unique()) if ID not in list(matched_rID)]
    added_gdf = setgdf[setgdf.Record_ID.isin(addedID)]
    if export:
        if missed_gdf.shape[0] > 0:
            selcols = ['wdID', 'trsqq', 'parcel_id', 'notes', 'lots', 'lot', 'ORTaxlot', 'record_ID', 'geometry']
            missed_gdf[~missed_gdf.geometry.isnull()][selcols].to_file(
                os.path.join(
                    INPATH,
                    f'\\{outfolder}\\',
                    f'missed_records_in_{setID}_res.shp'), 
                driver='ESRI Shapefile')
        if added_gdf.shape[0] > 0:
            added_gdf.rename(columns={'wetdet_delin_number': 'wdID', 
                      'address_location_desc':'loc_desc', 
                      'Coord_Source': 'coord_src',
                      'DocumentName':'doc_name',
                      'DecisionLink':'doc_link',
                      'is_batch_file':'isbatfile',
                      'status_name': 'status_nm',
                      'received_date':'receiveddt', 
                      'response_date':'responsedt',
                      'reissuance_response_date':'reissuance',
                      'Match_found':'matchfound', 
                      'Manual_note': 'notes',
                      'Edits_Complete': 'edits', 
                      'Shape_Length':'Shp_Length'           
                      }, inplace=True)
            added_gdf.to_file(
                os.path.join(
                    INPATH,
                    f'\\{outfolder}\\',
                    f'added_records_in_{setID}_res.shp'), 
                driver='ESRI Shapefile')
        
    return missed_match_ID, missed_gdf, addedID, added_gdf, missed_ID

def records_with_lots(gdf, setID, nm_to_add, c='Y'):
    """
    get wd records with lots
    input gdf from match_wd_data_with_taxlot
    return the unmatched records with/without taxlot IDs
    """
    unmatched_wd_df = report_unmatched(gdf, setID, nm_to_add)
    double_check = unmatched_wd_df[unmatched_wd_df.missinglot == c]
    return double_check

def review_with_lots(df, setID, all_taxlot, nm_to_add):
    """
    review the unmatched records with taxlot IDs to rematch with corrected data info
    input df is the reindexed dataframe from combined_reindexed_data
    return the matched taxlot, data frame without lots, the record IDs to map
    """
    n_gdf = match_wd_data_with_taxlot(df, setID, all_taxlot)
    wd_df = combine_wd_tables(setID, nm_to_add)
    df_wo_lots = records_with_lots(n_gdf, setID, nm_to_add, c='Y')
    wo_lots_ID = df_wo_lots.record_ID.unique()
    to_map_rID = [rID for rID in wd_df.record_ID.values if rID not in n_gdf.record_ID.unique()]
    unmatched_wlots_ID = [rID for rID in to_map_rID if rID not in wo_lots_ID]
    wlots_df = wd_df[wd_df.record_ID.isin(unmatched_wlots_ID)]  
    return n_gdf, wlots_df, to_map_rID

def review_mapped_data(df, setID, all_taxlot, nm_to_add, export=False):
    """
    get the digitized records from previous work
    input df is the reindexed dataframe from combined_reindexed_data
    return mapped data (and missed matched data) from the existing
    """
    setgdf = gpd.read_file(
        os.path.join(INPATH, 'GIS', 'Join_Statewide.gdb'),
        layer=f'WD_{setID}_Combined')
    to_map_rID = review_with_lots(df, setID, all_taxlot, nm_to_add)[2]
    mapped = setgdf[setgdf.Record_ID.isin(to_map_rID)]
    if export:
        mapped = mapped.rename(columns={'wetdet_delin_number': 'wdID', 
                      'address_location_desc':'loc_desc', 
                      'Coord-Source': 'coord_src',
                      'DocumentName':'doc_name',
                      'DecisionLink':'doc_link',
                      'is_batch_file':'isbatfile',
                      'status_name': 'status_nm',
                      'received_date':'receiveddt', 
                      'response_date':'responsedt',
                      'reissuance_response_date':'reissuance', 
                      'Match_found':'matchfound', 
                      'Manual_note': 'notes',
                      'Edits_Complete': 'edits'
                      })
        selcols = list(mapped.columns[list(map(lambda x: x <= 10, list(map(len, list(mapped.columns)))))])
        mapped[selcols].to_file(
            os.path.join(
                INPATH,
                f'\\{outfolder}\\',
                f'mapped_in_{setID}.shp'), 
            driver='ESRI Shapefile')
    return mapped

def check_corrected_data(df, setID, all_taxlot, nm_to_add, export=False):
    """
    update match on the corrected data
    input df is reindexed from combined_reindexed_data
    return matched data with corrected info
    """ 
    corrected = pd.read_excel(
        os.path.join(
            INPATH,
            'DSL data originals', 'Data corrections feedback to DSL', 
            f'DSL Database corrections {setID}.xlsx'))
    corrected.rename(columns={'trsqq': 'cor_trsqq',
                              'parcel_id':'cor_parcel_id'}, inplace=True)
    res = review_with_lots(df, setID, all_taxlot, nm_to_add)
    gdf = res[0]
    df_wlots = res[1]
    IDs1 = check_duplicates(df_wlots.wetdet_delin_number.values)[0]
    IDs2 = check_duplicates(corrected.wetdet_delin_number.values)[0]
    comIDs = [ID for ID in IDs1 if ID in IDs2]
    wdID_to_check = [wdID for wdID in df_wlots.wetdet_delin_number.values if wdID not in corrected.wetdet_delin_number.values]
    cols1 = df_wlots.columns
    cols2 = corrected.columns
    comcols = [col for col in cols1 if col in cols2]
    cols = [col for col in df_wlots.columns if col not in comcols]
    cols.append('wetdet_delin_number')
    df_wlots = df_wlots[cols]
    if len(comIDs) > 0:
        df_wlots = df_wlots[~df_wlots.wetdet_delin_number.isin(comIDs)]
        corrected = corrected[~corrected.wetdet_delin_number.isin(comIDs)]
    cor_df = corrected.merge(df_wlots, on = 'wetdet_delin_number')
    cor_df = cor_df.drop(['trsqq', 'parcel_id'], axis=1)
    cor_df.rename(columns={'cor_trsqq': 'trsqq', 
                          'cor_parcel_id': 'parcel_id',
                          'year': 'Year'}, inplace=True)
    df_wlots_to_check = df_wlots[df_wlots.wetdet_delin_number.isin(wdID_to_check+comIDs)]
    cor_df.loc[:, 'county'] = cor_df.county.apply(lambda x: string.capwords(x))
    cor_df = reindex_data(cor_df)
    cor_df_re = match_wd_data_with_taxlot(df=cor_df, setID=setID, all_taxlot=all_taxlot)
    collist = list(gdf.columns)
    collist.remove('recordID')
    ndf = gdf.append(cor_df_re[collist])
    if export:
        ngdf = gpd.GeoDataFrame(ndf, crs="EPSG:2992", geometry='geometry')
        selcols = ['wdID', 'trsqq', 'parcel_id', 'notes', 'lots', 'lot', 'ORTaxlot', 'record_ID', 'geometry']
        ngdf[~ngdf.geometry.isnull()][selcols].to_file(
            os.path.join(
                INPATH,
                '\\output\\matched\\',
                f'matched_records_{setID}.shp'),
            driver='ESRI Shapefile')
    return ndf, cor_df, comIDs, df_wlots_to_check

Tdict = dict(zip(['.25S', '.50S', '.50N', '.75S'], ['Z', 'V', 'Y', 'X']))
Rdict = dict(zip(['.25E', '.50E', '.50W', '.75E'], ['Z', 'V', 'Y', 'X']))

def taxlot2trsqq(x):
    """
    convert taxlot to trsqq
    """ 
    #print(x)
    rcs = ['.25', '.50', '.75'] # rcs = ratio codes
    if any([c in x for c in rcs]):
        if (x[4:7] in rcs) and (x[10:13] in rcs):
            if x[4:7] == x[10:13]:
                x1 = x.split(x[4:7])
            else:
                x0 = x.split(x[4:7])
                x1 = x0[1].split(x[10:13])
                x1.insert(0, x0[0]) 
            tcode = x[4:7] + x1[1][0]
            rcode = x[10:13] + x1[2][0]
            p1 = x1[0][2:4] + Tdict[tcode] + x1[1] + Rdict[rcode] + x1[2][0]
            if x1[2][1:5] == '0000':
                res = p1
            elif x1[2][3:5] == '00':
                res = p1 + x1[2][1:3]
            elif x1[2][4] == '0':
                res = p1 + x1[2][1:4]
            else:
                res = p1 + x1[2][1:5]

        elif x[4:7] in rcs:
            x1 = x.split(x[4:7])
            x2 = x1[1].split('.00')
            tcode = x[4:7] + x2[0][0]
            p2 = x1[0][2:4] + Tdict[tcode] + x2[0] + x2[1][0]
            if x2[1][1:5] == '0000':
                res = p2
            elif x2[1][3:5] == '00':
                res = p2 + x2[1][1:3]
            elif x2[1][4] == '0':
                res = p2 + x2[1][1:4]
            else:
                res = p2 + x2[1][1:5]
           
        elif x[10:13] in rcs:
            x1 = x.split(x[10:13])
            x2 = x1[0].split('.00')
            rcode = x[10:13] + x1[1][0]
            p3 = x2[0][2:4] + x2[1] + Rdict[rcode] + x1[1][0]
            if x1[1][1:5] == '0000':
                res = p3
            elif x1[1][3:5] == '00':
                res = p3 + x1[1][1:3]
            elif x1[1][4] == '0':
                res = p3 + x1[1][1:4]
            else:
                res = p3 + x1[1][1:5]
    elif any([s in x for s in ['.0S', '.0W', '.0E', '.0N']]):
        res = x[2:4] + x[6:9] + x[11:16]
                
    else:
        p4 = x[2:4] + x[7:10] + x[13]
        if x[14:18] == '0000':
            res = p4
        elif x[16:18] == '00':
            res = p4 + x[14:16]
        elif x[17] == '0':
            res = p4 + x[14:17]
        else:
            res = p4 + x[14:18]
            
    return res
    
# update match with similar trsqq
def get_trsqq_list():
    """
    get trsqq list, dictionary, and dataframe
    """    
    taxlotIDs_to_search = ALL_TXID
    taxlotIDs_cleaned = remove_duplicates([txID for txID in taxlotIDs_to_search if len(txID) == 29])
    trsqq = list(map(lambda x: taxlot2trsqq(x), taxlotIDs_cleaned))
    trsqq_dict = dict(zip(trsqq, taxlotIDs_cleaned))
    df = pd.DataFrame({'ORTaxlot':taxlotIDs_cleaned, 'trsqq':trsqq})
    with open(os.path.join(INPATH, "trsqq_list.pkl"), "wb") as f:
        pickle.dump(trsqq, f)
    with open(os.path.join(INPATH, "trsqq_dict.pkl"), "wb") as f:
        pickle.dump(trsqq_dict, f)
    df.to_csv(os.path.join(INPATH, "trsqq_df.csv"), index=False)
    return trsqq, trsqq_dict, df

def get_maybe_taxlot(trsqq_to_check):
    """
    get the maybe taxlot from the trsqq_to_check
    this is a test function
    """    
    closematch = difflib.get_close_matches(trsqq_to_check, TRSQQ)
    trsqq_matched = remove_duplicates(closematch)
    checktaxlot = [*map(TRSQQ_DICT.get, trsqq_matched)]
    string_to_search = trsqq_matched[0][0:8]
    search_res = remove_duplicates([i for i in remove_duplicates(TTDF.trsqq) if string_to_search in i])
    if len(checktaxlot) == 1 and len(search_res)==1:
        values = [i for i in TRSQQ_DICT if TRSQQ_DICT[i]==checktaxlot[0]]
    else:
        values = search_res
    return values

def reorganize_tocheck(tocheck_df):
    """ 
    reorganize the unmatched wd records to search for a potential match
    input is the output of check_corrected_data: df_wlots_to_check
    trsqq_n: new trsqq
    n_trsqq: number of possibly corrected trsqq values from the search
    this is a test function
    """
    tocheck_df.loc[:, 'trsqq_n'] = tocheck_df.loc[:, 'trsqq'].apply(lambda x: get_maybe_taxlot(x))
    tocheck_df.loc[:, 'n_trsqq'] = tocheck_df.loc[:, 'trsqq_n'].apply(lambda x: len(x))
    torematch_df = tocheck_df[tocheck_df.n_trsqq == 1]
    res = get_trsqq_list()
    trsqq_dict = res[1]
    torematch_df.loc[:, 'ORTaxlot'] = torematch_df.trsqq_n.apply(lambda x: [*map(trsqq_dict.get, x)][0])
    return tocheck_df, torematch_df

def run_Tier1(setID, nm_to_add, all_taxlot):
    """
    initial match
    """
    wd = combine_wd_tables(setID, nm_to_add)
    setdf = reindex_data(wd)
    # this might take a while
    start = time.time()
    # export = False
    setgdf = match_wd_data_with_taxlot(setdf, setID, all_taxlot, export=True)
    end = time.time()
    print(f'it took {round((end - start)/60, 0)} minutes to complete')
    unmatched_df = report_unmatched(setgdf, setID, nm_to_add, mute = False, export=True)
    return wd, setgdf, unmatched_df
