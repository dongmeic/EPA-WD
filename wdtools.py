import pandas as pd
import geopandas as gpd
import os
from os import walk
import zipfile
import re
import numpy as np
from itertools import chain
import time
import collections
import datetime
import string
import difflib
import pickle
import geopy
from geopy.geocoders import Nominatim
from urllib.request import urlopen
import io
import requests
from PyPDF2 import PdfReader, PdfWriter
import fiona
import webbrowser
import time
import googlemaps
import json
import openpyxl
from collections import Counter
from datetime import date
from win32com.client import Dispatch
from shapely.validation import make_valid
import warnings
from shapely.errors import ShapelyDeprecationWarning
from pyproj import Transformer
from random import sample
from shapely.geometry import shape
import glob
from PIL import Image
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
google_key=json.load(open('config/keys.json'))['google_maps']['APIKEY']

inpath = r'L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping'
outpath = inpath + '\\output'
wdpath = inpath + '\\DSL data originals'
txpath = inpath + '\\GIS\\ORMAP_data\\ORMAP_Taxlot_Years'
yearstart = 2016
yearend = 2023
#outfolder = 'test'
outfolder = 'output\\matched'
# create a spreadsheet to create a dictionary for the match between county name and code
cnt_ID = pd.read_excel(inpath+'\\notes\\CNT_Code.xlsx')
# create a dictionary to look up county code
cnt_dict = dict(zip(cnt_ID.COUNTY, cnt_ID.ID))
trsqq_correction_dict = dict(zip(list(range(0, 6)), ['township number', 'township direction', 'range number', 'range direction', 'section number', 'QQ']))
OR_counties = list(cnt_dict.keys())
nm2add = [0, 1420, 2143, 2878, 3932, 4370, 4972]
selectedvars = ['wetdet_delin_number', 'trsqq', 'parcel_id','address_location_desc', 
           'city', 'county', 'site_name', 'site_desc','latitude',
           'longitude', 'DocumentName', 'DecisionLink','is_batch_file',
           'status_name', 'received_date', 'response_date','reissuance_response_date', 
           'project_id', 'site_id', 'SetID','record_ID', 'ORMapNum']
varlist = selectedvars + ['geometry', 'code']
transformer = Transformer.from_crs("EPSG:2992", "EPSG:4326")
coldict = {'wetdet_delin_number': 'wdID', 
           'address_location_desc':'loc_desc', 
           'Coord-Source':'CordSource',
           'DocumentName':'doc_name',
           'DecisionLink':'doc_link',
           'is_batch_file':'isbatfile',
           'status_name': 'status_nm',
           'received_date':'receiveddt', 
           'response_date':'responsedt',
           'reissuance_response_date':'reissuance'}
issuepath = inpath + '\\GIS\\ArcGIS Pro Project\\DataReview\\issueIDs.gdb'

def read_trsqq():
    """
    read trsqq list, dictionary, and dataframe
    """    
    with open(os.path.join(inpath, "trsqq_list.pkl"), "rb") as f:
        trsqq = pickle.load(f)
    with open(os.path.join(inpath, "trsqq_dict.pkl"), "rb") as f:
        trsqq_dict = pickle.load(f) 
    df = pd.read_csv(os.path.join(inpath, "trsqq_df.csv"))
    return trsqq, trsqq_dict, df

trsqq, trsqq_dict, ttdf = read_trsqq()
tid_dst = [tid for tid in ttdf.ORTaxlot.unique() if any(substring in tid for substring in ['--D', '--S', '--T'])]
tid_dst_0 = list(map(lambda x: re.split("--", x)[0], tid_dst))
tid_dst_1 = list(map(lambda x: re.split("--", x)[1], tid_dst))
tsq_dst = ttdf[ttdf.ORTaxlot.isin(tid_dst)].trsqq.unique()
cnts = gpd.read_file(inpath + "\\GIS\\Oregon_Counties.shp")

pdf_outpath = r'L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping\output\pdf'
with open(os.path.join(inpath, "ORTaxlot.pkl"), "rb") as f:
    all_txid = pickle.load(f)
with open(os.path.join(inpath, "ORMapIndex.pkl"), "rb") as f:
    all_mpidx = pickle.load(f)
with open(os.path.join(inpath, "ParticipCnt.pkl"), "rb") as f:
    cntlst = pickle.load(f)
pd.options.mode.chained_assignment = None

############################################### QAQC updates ###################################################

def update_QAQC_data(setID, wd, totcol, qaqc_col, export=True):
#     file = outpath + f'\\to_review\\re_mapping_{setID}.txt'
#     with open(file) as f:
#         remapIDs = f.readlines()
#     partialIDs = list(wd[~wd.notes.isnull()].wetdet_delin_number.unique())
#     unmatched = pd.read_csv(outpath+f'\\to_review\\unmatched_df_{setID}_2.csv')
#     unmatchedIDs = list(unmatched.wetdet_delin_number.unique()) 
#     issues = pd.read_excel(outpath + f'\\to_review\\{setID}_Mapping_Issues.xlsx')
#     issueIDs = list(issues.wetdet_delin_number.unique())
#     qaqcIDs = remapIDs[0].split(', ') + partialIDs + unmatchedIDs + issueIDs
    mapped = gpd.read_file(f'{outpath}//final//mapped_wd_{setID}.shp')
    qaqcIDs = wd[~wd.wetdet_delin_number.isin(mapped.wdID.unique())].wetdet_delin_number.unique()
    # count QAQC 
    qaqc_df = wd[wd.wetdet_delin_number.isin(qaqcIDs)][['county', 'wetdet_delin_number']].groupby(['county']).agg(lambda x: x.nunique()).reset_index().rename(columns={'wetdet_delin_number':'QAQC_count'})
    # total count
    total_df = wd[['county', 'wetdet_delin_number']].groupby(['county']).agg(lambda x: x.nunique()).reset_index().rename(columns={'wetdet_delin_number':'total_count'})
    # update the output
    gdf = gpd.read_file(inpath+'\\reporting\\WD_Counts.shp')
    sel = qaqc_df.county.unique()
    wdcnt_dict = total_df.set_index('county').to_dict(orient='dict')['total_count']
    qccnt_dict = qaqc_df.set_index('county').to_dict(orient='dict')['QAQC_count']
    gdf.loc[gdf.County.isin(sel), totcol] = gdf[gdf.County.isin(sel)].County.map(wdcnt_dict)
    gdf.loc[gdf.County.isin(sel), qaqc_col] = gdf[gdf.County.isin(sel)].County.map(qccnt_dict)
    gdf.loc[:, 'Tot_ToMap'] =  gdf.P1_ToMap + gdf.P2_ToMap + gdf.P3_ToMap
    gdf.loc[:, 'Tot_Count'] =  gdf.P1_2017_22 + gdf.P2_2008_16 + gdf.P3_1990_07
    gdf.loc[:, 'QAQC_Ratio'] = gdf.Tot_ToMap / gdf.Tot_Count
    gdf.loc[:, 'Auto_Comp'] = gdf.Tot_Count - gdf.Tot_ToMap
    if export:
        gdf.to_file(inpath+'\\reporting\\WD_Counts.shp')
    return gdf

def file_compress(inp_file_names, out_zip_file):
    """
    function : file_compress
    args : inp_file_names : list of filenames to be zipped
    out_zip_file : output zip file
    return : none
    assumption : Input file paths and this code is in same directory.
    """
    # Select the compression mode ZIP_DEFLATED for compression
    # or zipfile.ZIP_STORED to just store the file
    compression = zipfile.ZIP_DEFLATED
    print(f" *** Input File name passed for zipping - {inp_file_names}")

    # create the zip file first parameter path/name, second mode
    print(f' *** out_zip_file is - {out_zip_file}')
    zf = zipfile.ZipFile(out_zip_file, mode="w")
    
    try:
        for file_to_write in inp_file_names:
            # Add file to the zip file
            # first parameter file to zip, second filename in zip
            print(f' *** Processing file {file_to_write}')
            zf.write(file_to_write, file_to_write, compress_type=compression)

    except FileNotFoundError as e:
        print(f' *** Exception occurred during zip process - {e}')
    finally:
        # Don't forget to close the file!
        zf.close()

############################################### Taxlot Review ###################################################
def removeCountyNm(x, toRemove):
    """
    remove county name from strings
    """      
    tf = [tr in x for tr in toRemove]
    tr = [tr for tr in toRemove if tr in x]
    if any(tf):
        x=re.sub(tr[0], '', x)
    return x

def readGeoData(layer_file):
    """
    read geodata with a geometry check
    """ 
    try:
        gdf = gpd.read_file(layer_file)
    except ValueError:
        collection = list(fiona.open(layer_file,'r'))
        df1 = pd.DataFrame(collection)

        #Check Geometry
        def isvalid(geom):
            try:
                shape(geom)
                return 1
            except:
                return 0
        df1['isvalid'] = df1['geometry'].apply(lambda x: isvalid(x))
        df1 = df1[df1['isvalid'] == 1]
        collection = json.loads(df1.to_json(orient='records'))

        #Convert to geodataframe
        gdf = gpd.GeoDataFrame.from_features(collection)
#     if gdf.crs is None:
#         gdf.crs = "epsg:2992"
    return gdf

def readTaxlots(year, ORTxt_Only=False, export=False):
    """
    read taxlots prior to 2017 
    """
    frames = []
    if ORTxt_Only:
        listcols = ['ORTaxlot', 'geometry']
    else:
        listcols = ['County', 'Town', 'TownPart', 'TownDir','Range','RangePart','RangeDir', 'SecNumber', 'Qtr', 'QtrQtr', 'Anomaly', 'MapSufType', 'MapNumber', 'ORMapNum', 'Taxlot','MapTaxlot', 'ORTaxlot', 'geometry']
    selcols = list(map(lambda x: x.capitalize(), listcols))
    gdb = inpath+f'\\GIS\\ORMAP_data\\raw\\Taxlots{year}.gdb'
    if year in range(2016, 2018):      
        txlot = gpd.read_file(gdb, layer='Taxlots')
        if all([col in txlot.columns for col in listcols]):
            txlot = txlot[listcols]
        else:
            print(f"check the column names of the {year} taxlots!")
        frames.append(txlot)
    elif year in [2011, 2014, 2015]:
        lyrlist = fiona.listlayers(gdb)
        lyrsel = [x for x in lyrlist if x not in list(filter(lambda x : re.search(r"_prop_tbl|_TaxCode", x), lyrlist))]
        if year != 2011:
            lyrsel = [lyr for lyr in lyrsel if lyr in cntlst]
        else:
            lyrsel = [lyr for lyr in lyrsel if (lyr in cntlst) and (lyr != 'Jefferson')]
        toRemove = list(map(lambda x: x+'_', lyrsel))
        for lyr in sorted(lyrsel):
            print(lyr)
            txlot = gpd.read_file(gdb, layer=lyr)
            txlot = txlot.to_crs(epsg=2992)
            lst = list(map(lambda x: re.sub(r'[0-9]', '', x), txlot.columns))
            colnms = list(map(lambda x:removeCountyNm(x, toRemove), lst))
            txlot.columns = colnms
            if year != 2011:
                txt = 'taxlot__|taxlot_|Taxlots_|TAXLOT_|Taxlot_|taxlots_|TaxLot_'
                colsel = list(filter(lambda x: re.search(txt, x, re.IGNORECASE), colnms))
                colsel = unique([col for col in colsel if col not in ['Taxlots_F', 'Taxlots_Map_Taxlo', 'Taxlots_Accnum']]) + ['geometry']
                ncolnms = list(map(lambda x: re.sub(txt, '', x, re.IGNORECASE), colsel))
                txlot = txlot[colsel]
                txlot.columns = list(map(lambda x: x.capitalize(), ncolnms))
            else:
                if all([colnm in colnms for colnm in ['MapTaxlot', 'Maptaxlot']]):
                    colsel = [col for col in colnms if col not in check_duplicates(colnms)[0]+['Maptaxlot']]
                else:
                    colsel = [col for col in colnms if col not in check_duplicates(colnms)[0]]
                txlot = txlot[colsel]
                txlot.columns = list(map(lambda x: x.capitalize(), colsel))
            txlot = txlot[selcols]  
            txlot.columns = listcols
            frames.append(txlot)
    elif year == 2012:
        dir_list = os.listdir(inpath+f'\\GIS\\ORMAP_data\\raw\\Taxlots{year}')
        # Umatilla and Lane don't have the complete required information
        lyrs = [lyr for lyr in dir_list if lyr not in ['Umatilla']]
        for lyr in lyrs:
            print(lyr)
            path = inpath+f'\\GIS\\ORMAP_data\\raw\\Taxlots{year}\\{lyr}'
            txt = 'taxlot|Taxlot'
            files = os.listdir(path)
            filelst = list(filter(lambda x: re.search(txt, x, re.IGNORECASE), files))
            filesel = [file for file in filelst if 'Taxlots' not in file]
            file = filesel[0].split('.')[0]
            txlot = readGeoData(path+f'\\{file}.shp')
            if lyr == 'Lane':
                txlot['TownDir'] = txlot.ORTaxlot.apply(lambda x: separate_numbers_letters(str(x))[1][0])
            if txlot.crs is None:
                if lyr in ['Baker', 'Gilliam', 'Hood River', 'Jefferson', 'Columbia', 'Lincoln', 'Marion', 'Polk', 'Tilamook', 
                          'Union', 'Wasco', 'Jackson', 'Wheeler', 'Yamhill']:
                    txlot = txlot.set_crs('epsg:2913')
                    txlot = txlot.to_crs('epsg:2992')
                elif lyr in ['crook', 'curry', 'Josephine', 'Klamath', 'Lake', 'Lane']:
                    txlot = txlot.set_crs('epsg:2914')
                    txlot = txlot.to_crs('epsg:2992')
                elif lyr in ['Benton', 'Clackamas', 'Clatsop', 'Linn']:
                    txlot = txlot.set_crs('epsg:2994')
                    txlot = txlot.to_crs('epsg:2992')
                elif lyr in ['Deschutes', 'Harney', 'Jackson']:
                    txlot = txlot.set_crs('epsg:2270')
                    txlot = txlot.to_crs('epsg:2992')
                elif lyr == 'Umatilla':
                    txlot = txlot.set_crs('epsg:32126')
                    txlot = txlot.to_crs('epsg:2992')
                elif lyr in ['Multnomah', 'Washington']:
                    txlot = txlot.set_crs('epsg:102326')
                    txlot = txlot.to_crs('epsg:2992')
                else:
                    print(f'check the projected coordinate system for {lyr} in 2012')
            else:
                txlot = txlot.to_crs('epsg:2992')
            if all([col in txlot.columns for col in listcols]):
                txlot = txlot[listcols]
            elif all([col in list(map(lambda x: x.capitalize(), txlot.columns)) for col in selcols]):
                txlot.columns = list(map(lambda x: x.capitalize(), txlot.columns))
                txlot = txlot[selcols]  
                txlot.columns = listcols
            else:
                print(f"check the column names of the {year} taxlots in {lyr}!")
            frames.append(txlot)        
    elif year == 2009:
        layers = glob.glob(inpath+f'\\GIS\\ORMAP_data\\raw\\Taxlots{year}\\' + '*.shp')
        for layer in layers:   
            c = layer.split('Taxlots2009\\')[1].replace('_tax_09.shp' , '').replace('_Tax_09.shp' , '')
            cl = [cnt for cnt in cnt_ID.COUNTY.values if (cnt[0:4] == c) or (cnt[0:3] == c)]
            if cl[0] != 'Wasco': 
                #print(layer)
                gdf = gpd.read_file(layer)
                gdf = gdf.to_crs(epsg=2992)
                gdf.loc[:, 'County'] = cl[0]
                #print(gdf.columns)  
                if 'ORTaxlot' in gdf.columns:
                    gdf = gdf[['ORTaxlot', 'County', 'geometry']]
                    frames.append(gdf)
                else:
                    if any([col in gdf.columns for col in ['FIRST_ORTa', 'FIRST_ORTA', 'FIRST_ORMa']]):

                        gdf.rename(columns={'FIRST_ORTa': 'ORTaxlot',
                                            'FIRST_ORTA': 'ORTaxlot',
                                            'FIRST_ORMa': 'ORTaxlot'}, 
                                   inplace=True)
                    else:
                        gdf.rename(columns={'SIMAPTAX': 'ORTaxlot',
                                           'ORTAXLOT': 'ORTaxlot'}, 
                                   inplace=True)
                    if 'ORTaxlot' in gdf.columns:
                        gdf = gdf[['ORTaxlot', 'County', 'geometry']]
                        frames.append(gdf)
    else:
        print("check the year!")
    gdf = pd.concat(frames, ignore_index=True)
    gdf.loc[:, 'Year'] = year
    if export:
        gdf = gdf[~gdf.geometry.isnull()]
        if isinstance(gdf, pd.DataFrame):
            gdf = gpd.GeoDataFrame(gdf, crs="EPSG:2992", geometry="geometry")
        gdf.to_file(inpath + f'\\GIS\\ORMAP_data\\2009_2015\\ORTaxlots{year}.shp')
    return gdf

def readMapIndex(year):
    colnms = ['County', 'ORMapNum', 'geometry']
    if year == 2012:
        frames = []
        dir_list = os.listdir(txpath + f'\\Taxlots{year}')
        for lyr in dir_list:
            #print(lyr)
            path = txpath + f'\\Taxlots{year}\\{lyr}'
            files = os.listdir(path)
            filelst = list(filter(lambda x: re.search('mapindex|MIndex', x, re.IGNORECASE), files))
            file = filelst[0].split('.')[0]
            layer_file = path+f'\\{file}.shp'
            gdf = readGeoData(layer_file)
            if all([col in gdf.columns for col in colnms]):
                gdf = gdf[colnms]
            else:
                gdf.columns = list(map(lambda x: x.capitalize(), gdf.columns))
                if all([col in gdf.columns for col in ['County', 'Ormapnum', 'Geometry']]):
                    gdf.rename(columns={'Ormapnum': 'ORMapNum', 'Geometry': 'geometry'}, inplace=True)
                    gdf = gdf[colnms]
                elif 'County' not in gdf.columns:
                    gdf.loc[:, 'County'] = int(cnt_dict[lyr])
                    if 'First_orma' in gdf.columns:
                        gdf.rename(columns={'First_orma': 'ORMapNum', 'Geometry': 'geometry'}, inplace=True)
                    else:
                        gdf.rename(columns={'Ormapnum': 'ORMapNum', 'Geometry': 'geometry'}, inplace=True)
                    gdf = gdf[colnms]
                else:
                    print(f"check the column names of the {year} taxmap in {lyr}!")
            frames.append(gdf)
        gdf = pd.concat(frames, ignore_index=True)
    elif year == 2011:
        gdb = txpath+f'\\Taxlots{year}.gdb'
        gdf = gpd.read_file(gdb, layer='mapIndex')
        gdf = gdf[colnms]
    else:
        print("check the year!")
    gdf.loc[:, 'Year'] = year
    return gdf

################################################ Deliverable ####################################################

def format_gdf_provided(gdf, wdID):
    """
    format gdf provided to be in the same format with mapped records
    """      
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
    """
    combine all the corrected DSL tables into one dataframe
    """  
    frames = []
    if setID:
        wd_df = pd.read_csv(wdpath+f'\\Corrected_by_Set\\Set{num}_2017-20220601.csv')
    else:
        for i in range(num):
            wd_dt = pd.read_csv(wdpath+f'\\Corrected_by_Set\\Set{i+1}_2017-20220601.csv')
            frames.append(wd_dt)
        wd_df = pd.concat(frames, ignore_index=True)
    return wd_df

def export_wd_gdf_by_record(gdf, outnm):
    """
    export geodataframe with the original wd data
    """      
    gdf = gdf[varlist]
    gdf['received_date'] = gdf['received_date'].dt.strftime("%Y-%m-%d")
    gdf['response_date'] = gdf['response_date'].dt.strftime("%Y-%m-%d")
    gdf['lat'], gdf['lon'] = transformer.transform(gdf.centroid.x, gdf.centroid.y)
    gdf = gdf.rename(columns=coldict)
    try:
        gdf.to_file(f'{outpath}\\test\\{outnm}.shp')
    except RuntimeError:
        gdf['geometry'] = gdf.geometry.buffer(0)
    gdf.to_file(f'{outpath}\\test\\{outnm}.shp')
    return gdf
        
def split_SA_by_rid_in_df(wd_df, sa_gdf_all, all_mapIdx, all_taxlot, em_wids, export=False, outnm='example_data', review=False):
    """
    split study area polygons where multiple record IDs exist;
    rid is record ID;
    em_wids is the example WD IDs;
    wd_df is the dataframe that includes the example WD IDs;
    sa_gdf_all is the geodataframe that includes the example WD IDs;
    all_mapIdx is the combined geodataframe of map index from all years;
    all_taxlot is the combined geodataframe of taxlots from all years;
    return the combined geodataframe from split_WD_to_records
    """
    wd_df_s = wd_df[wd_df.wetdet_delin_number.isin(em_wids)]
    wd_df_s['ORMapNum'] = wd_df_s[['county', 'trsqq']].apply(lambda row: create_ORMapNm(ct_nm=row.county, trsqq=row.trsqq), axis = 1)
    sa_gdf_s = sa_gdf_all[sa_gdf_all.wdID.isin(em_wids)]
    wdID_list, n=check_duplicates(wd_df_s.wetdet_delin_number.values)
    if n>0:
        frames = []
        for wid in wdID_list:
            if wid in sa_gdf_s.wdID.values:
                #print(wid)
                ORmn = wd_df_s[wd_df_s.wetdet_delin_number==wid].ORMapNum.values
                yr = wid[2:6]
                mapIdx = all_mapIdx[(all_mapIdx.ORMapNum.isin(ORmn)) & (all_mapIdx.year==yr)]
                taxlot = all_taxlot[all_taxlot.year==yr]
                if review:
                    out = split_WD_to_records(df=wd_df_s, gdf=sa_gdf_s, wdID=wid, mapindex=mapIdx, taxlots=taxlot, review=True)
                else:
                    out = split_WD_to_records(df=wd_df_s, gdf=sa_gdf_s, wdID=wid, mapindex=mapIdx, taxlots=taxlot)
                frames.append(out)
        df = pd.concat(frames, ignore_index=True)
        gdf1 = gpd.GeoDataFrame(df, crs="EPSG:2992", geometry='geometry')
        wdID_sdf = wd_df_s[~wd_df_s.wetdet_delin_number.isin(wdID_list)]
        wdIDList = wdID_sdf.wetdet_delin_number.values
        sa_sgdf = sa_gdf_s[sa_gdf_s.wdID.isin(wdIDList)]
        sa_sgdf.rename(columns={'wdID':'wetdet_delin_number'}, inplace=True)
        gdf2 = wdID_sdf.merge(sa_sgdf[['code', 'wetdet_delin_number', 'geometry']], on='wetdet_delin_number')
        gdf = pd.concat([gdf1[varlist], gdf2[varlist]])
    else:
        sa_gdf_s.rename(columns={'wdID':'wetdet_delin_number'}, inplace=True)
        gdf = wd_df_s.merge(sa_gdf_s[['code', 'wetdet_delin_number', 'geometry']], on='wetdet_delin_number')
        gdf = gdf[varlist]
        gdf = gpd.GeoDataFrame(gdf, crs="EPSG:2992", geometry='geometry')
    #gdf['received_date'] = gdf['received_date'].dt.strftime("%Y-%m-%d")
    #gdf['response_date'] = gdf['response_date'].dt.strftime("%Y-%m-%d")
    gdf['lat'], gdf['lon'] = transformer.transform(gdf.representative_point().x, gdf.representative_point().y)
    if export:
        gdf = gdf.rename(columns=coldict)
        try:
            gdf.to_file(f'{outpath}\\test\\{outnm}.shp')
        except RuntimeError:
            gdf['geometry'] = gdf.geometry.buffer(0)
            gdf.to_file(f'{outpath}\\test\\{outnm}.shp')    
    return gdf
           
def read_all_mapIdx(exportID=False):
    """
    combine mapIndex from all years
    """
    frames = []
    for year in range(yearstart, yearend):
        mapIdx_dt = read_mapIdx(year)
        mapIdx_dt['year'] = str(year)
        frames.append(mapIdx_dt[['year', 'ORMapNum', 'geometry']])
    df = pd.concat(frames, ignore_index=True)
    gdf = gpd.GeoDataFrame(df, crs="EPSG:2992", geometry='geometry')
    if exportID:
        with open(os.path.join(inpath, "ORMapIndex.pkl"), "wb") as f:
            pickle.dump(list(gdf.ORMapNum.unique()), f) 
    return gdf
    
def read_mapIdx(year):
    """
    read mapIndex from one year
    """
    mapindex = gpd.read_file(inpath+f'\\GIS\\ORMAP_data\\ORMAP_Taxlot_Years\\Taxlots{year}.gdb', 
                             layer='MapIndex')
    return mapindex

def replace_geometry(gdf):
    """
    to replace the geometry of SA ploygons with manual review
    gdf is the geodataframe that contains the SA polygons to replace
    rid is record ID
    """
    revpath = inpath + '\\GIS\\ArcGIS Pro Project\\DataReview\\DataReview.gdb'
    newplys = [lyr for lyr in fiona.listlayers(revpath) if 'rid' in lyr]
    frames = []
    for nply in newplys:
        newply = gpd.read_file(revpath, layer=nply)
        rid = int(nply.replace('rid', ''))
        newply['record_ID'] = rid
        frames.append(newply[['record_ID', 'geometry']])
    sa_df = pd.concat(frames, ignore_index=True)
    sa_gdf = gpd.GeoDataFrame(sa_df, geometry='geometry')
    selrids = sa_df.record_ID.unique()
    df = gdf.drop(columns=['geometry'])
    gdf1 = df[df.record_ID.isin(selrids)].merge(sa_gdf,on='record_ID')
    gdf2 = pd.concat([gdf[~gdf.record_ID.isin(selrids)], gdf1], ignore_index=True)
    return gdf2
    
def split_WD_to_records(df, gdf, wdID, mapindex, taxlots, review=False):
    """
    this function splits the WD SA ploygons to ploygons by records
    df is the dataframe that contains the selected WD ID and wetdet_delin_number in the columns
    gdf is the geodataframe that contains the selected WD ID
    mapindex is the taxmap geodataframe of the year
    taxlots are the taxlots of the year
    """
    ndf = df[df.wetdet_delin_number==wdID]
    cnts = ndf.county.unique()
    if (len(cnts) > 1) and (review==False):
        print(f"{wdID} crosses counties!")
        return None
    elif cnts.any() not in OR_counties:
        print(f"Check the county name {cnts[0]}!")
        return None
    else:
        setID = ndf.SetID.values[0]
        #print(f"{wdID} is in County {cnts[0]} in Set {setID}...")
        if 'ORMapNum' not in ndf.columns:
            ndf['ORMapNum'] = ndf[['county', 'trsqq']].apply(lambda row: create_ORMapNm(ct_nm=row.county, trsqq=row.trsqq), axis = 1)
        trsqq_list, n = check_duplicates(ndf.trsqq.values)
        ngdf = split_WD_to_taxmaps(df=ndf, gdf=gdf, wdID=wdID, mapindex=mapindex)
        if n > 0:
            frames = []
            for trsqq in trsqq_list:
                inter = split_taxmap_to_records(df=ndf, gdf=ngdf, trsqq=trsqq, taxlots=taxlots)
                frames.append(inter)
            # idf is intersection data frame
            idf = pd.concat(frames, ignore_index=True)
            idf = idf[['geometry', 'record_ID']].merge(ndf[selectedvars],on='record_ID', how='left')
            idf['code'] = ngdf.code.values[0]
            # odf is the other dataframe that excludes intersection taxmaps
            ogdf=ngdf[~ngdf.ORMapNum.isin(idf.ORMapNum)][['ORMapNum', 'geometry', 'code']]
            odf=ndf[~ndf.ORMapNum.isin(idf.ORMapNum)][selectedvars]
            odf=ogdf.merge(odf, on='ORMapNum', how='left')
            out=pd.concat([idf[varlist], odf[varlist]])
        else:
            out=ndf.merge(ngdf[['ORMapNum', 'geometry', 'code']], on='ORMapNum', how='left')
            out=out[varlist]
        out = gpd.GeoDataFrame(out, geometry='geometry')
        out = out.dissolve('record_ID')
        out['record_ID'] = out.index
        out.reset_index(drop=True, inplace=True)
        return out

def split_WD_to_taxmaps(df, gdf, wdID, mapindex):
    """
    this function splits the WD SA ploygons by taxmap
    gdf is the geodataframe that contains the selected WD ID
    mapindex is the taxmap geodataframe of the year
    """
    selmid = df[df.wetdet_delin_number==wdID].ORMapNum.unique()
    gdf = gdf[gdf.wdID==wdID]
    selmapindex = mapindex[['ORMapNum','geometry']][mapindex.ORMapNum.isin(selmid)]
    try:
        inter = gpd.overlay(gdf, selmapindex, 
                    how='intersection', keep_geom_type=False)
    except NotImplementedError:
        gdf['geometry'] = gdf.geometry.buffer(0)
        inter = gpd.overlay(gdf, selmapindex, 
                    how='intersection', keep_geom_type=False)
    inter = inter.dissolve('ORMapNum')
    inter['ORMapNum'] = inter.index
    inter.reset_index(drop=True, inplace=True)
    return inter

def split_taxmap_to_records(df, gdf, trsqq, taxlots):
    """
    this function is applied when the same taxmap has multiple records
    df is the dataframe of the selected WD ID
    gdf is the geodataframe of the selected taxmaps, inter from split_WD_to_taxmaps
    trsqq is the selected trsqq that appears in multiple record IDs
    taxlots are the taxlots of the year
    """
    df = df[df.trsqq==trsqq]
    rdf = reindex_data(df)
    taxlots = taxlots[taxlots.ORTaxlot.isin(rdf.ORTaxlot.unique())]
    t_df = rdf[['ORTaxlot','record_ID']].merge(taxlots[['ORTaxlot', 'geometry']], 
                                                     on='ORTaxlot', 
                                                     how='left')
    t_gdf = gpd.GeoDataFrame(t_df, geometry='geometry')
    t_gdf = t_gdf.dissolve('record_ID')
    t_gdf['record_ID'] = t_gdf.index
    i_gdf = gdf[gdf.ORMapNum==df.ORMapNum.unique()[0]]
    i_gdf = i_gdf.dissolve('wdID')
    i_gdf['wdID'] = i_gdf.index
    i_gdf.reset_index(drop=True, inplace=True)
    inter = gpd.overlay(i_gdf, t_gdf, how='intersection', 
                        keep_geom_type=False)
    return inter
    
def create_ORMapNm(ct_nm, trsqq):
    """
    return ORMap number based on county name and trsqq
    """
    part1 = str(int(cnt_dict[ct_nm])).zfill(2) + convert_trsqq(trsqq)
    mpidx = part1 + '--0000'
    if mpidx not in all_mpidx:
        if part1 in tid_dst_0:
            for mid in [part1+f'--{x}000' for x in ['D', 'S', 'T']]:
                if mid in all_mpidx:
                    mpidx = mid     
    return mpidx
    
def get_all_wd(num, raw=False):
    """
    combine all the original DSL tables into one dataframe
    """    
    frames = []
    for i in range(num):
        if raw:
            wd_dt = combine_wd_tables(setID='Set00'+str(i+1), nm_to_add=nm2add[i], raw=True)
        else:
            wd_dt = combine_wd_tables(setID='Set00'+str(i+1), nm_to_add=nm2add[i], raw=False)
        wd_dt['SetID'] = i+1
        frames.append(wd_dt)
    wd_df = pd.concat(frames, ignore_index=True)
    return wd_df

def get_added_SA():
    path = inpath + "\\GIS\\ArcGIS Pro Project\\DataReview\\added.gdb"
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
        file1 = outpath + f'\\final\\mapped_wd_Set00{i+1}.shp'
        file2 = outpath + f'\\final\\Set00{i+1}_mapped_wd.shp'
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
    df['ORMapNum'] = df[['county', 'trsqq']].apply(lambda row: create_ORMapNm(ct_nm=row.county, trsqq=row.trsqq), axis = 1)
    for wid in wdlist:
        #print(wid)
        df_s = df[df.wetdet_delin_number==wid]
        gdf_s = gdf[gdf.wdID==wid]
        ORmn = df_s.ORMapNum.values
        yr = wid[2:6]
        mapIdx = mapindex[(mapindex.ORMapNum.isin(ORmn)) & (mapindex.year==yr)]
        exp_gdf = split_WD_to_taxmaps(df=df_s, gdf=gdf_s, wdID=wid, mapindex=mapIdx)
        exp_gdf['lat'], exp_gdf['lon'] = transformer.transform(exp_gdf.representative_point().x, exp_gdf.representative_point().y)
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
        sa_gdf=sa_gdf.rename(columns=coldict)
        try:
            sa_gdf.to_file(os.path.join(inpath, "output", "final", f"{outnm}.shp"), index=False)
        except RuntimeError:
            sa_gdf['geometry'] = sa_gdf.geometry.buffer(0)
            sa_gdf.to_file(os.path.join(inpath, "output", "final", f"{outnm}.shp"), index=False)
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
    with open(os.path.join(inpath, f"{setid}_mapped.pkl"), "rb") as f:
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
    dat = pd.read_excel(file)
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
    with open(os.path.join(inpath, f"{setID}_{lstnm}.pkl"), "wb") as f:
            pickle.dump(lst, f)
    
def review_loop_r1(setID=None, wdid_list=None, df=None, partial=False, idx=False, wd_id=None, wddf=None, plot=False, gdf=None):
    """
    loop through the unmatched records and check the original records
    """
    toadd=[]
    if not partial:
        df = pd.read_csv(os.path.join(inpath + f'\\output\\to_review\\unmatched_df_{setID}_r1_N.csv'))
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
        print(f'{round(((j/n)*100),1)}% digitized, {n-j} records remained, expected to be done in about {int(((n-j)*0.25)+0.5)} hours...')
        print(wdID)
        if idx:
            if df is not None:
                print(f'index = {df[df.wetdet_delin_number==wdID].index[0]+1}')
                print(check_unmatched_r1(wdID = wdID, df = df))
            else:
                print(f'index = {wdid_list.index(wdID)+1}')
                print(check_unmatched_r1(wdID = wdID, df = wddf))
        if plot:
            gdf[gdf.wdID == wdID].plot()
            plt.savefig('SA_plot.jpg')
            img = Image.open('SA_plot.jpg')
            img.show()
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
    df = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_N_0.csv'))
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
    revpath = inpath + f'\GIS\ArcGIS Pro Project\DataReview\Set{setID}.gdb'
    partial = gpd.read_file(revpath, layer=f'{setID}_partial')
    mapped1 = list(partial.wdID.unique())
    mapped0 = [lyr for lyr in fiona.listlayers(revpath) if (lyr not in [f'Set{setID}_wo_lot', f'Set{setID}_partial']) and ('L' not in lyr)]
    mapped2 = list(map(lambda x: x.replace('_', '-'), mapped0))
    mapped = mapped1 + mapped2
    matched = gpd.read_file(inpath + f'\\output\matched\matched_records_Set{setID}_edited.shp')
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
    revpath = inpath + f'\GIS\ArcGIS Pro Project\DataReview\{setID}.gdb'
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
        revpath = inpath + f'\GIS\ArcGIS Pro Project\DataReview\{setID}.gdb'
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

def read_text_file(file):
    """
    read text file to get a list of WD IDs
    """
    if os.path.exists(file):
        with open(file) as f:
            IDs = f.readlines()
            if len(IDs) == 1:
                IDs = IDs[0].split(', ')
    else:
        IDs = []
    return IDs

def combine_matched_digitized(setID, editedIDs, nm_to_add, export=True):
    """
    Combine the edited matched records, digitized partial taxlots, taxlots without lot IDs, and the list of issue IDs
    """
    revpath = inpath + f'\GIS\ArcGIS Pro Project\DataReview\{setID}.gdb'
    # get separated feature files
    mapped0 = [lyr for lyr in fiona.listlayers(revpath) if (lyr not in [f'{setID}_wo_lot', f'{setID}_partial']) and ('L' not in lyr)]
    matched = gpd.read_file(inpath + f'\\output\matched\matched_records_{setID}_edited.shp')
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
    issues = pd.read_csv(os.path.join(inpath, "output", "to_review", f"{setID}_Mapping_Issues.csv"))
    # exclude the ones that have issues
    issueIDs = list(issues.wetdet_delin_number.unique())
    reviewedIDs = read_text_file(file = outpath+f"\\matched\\{setID}_reviewed.txt")
    issueIDs = [iID for iID in issueIDs if iID not in reviewedIDs]
    withdrawnIDs = read_text_file(file = outpath+f"\\matched\\{setID}_not_mapped.txt")
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
    wd = combine_wd_tables(setID=setID, nm_to_add=nm_to_add)
    unmatchedIDs = [wdID for wdID in wd.wetdet_delin_number.unique() if wdID not in final_gdf.wdID.unique()]
    toCheck = [ID for ID in unmatchedIDs if ID not in issueIDs]
    digitized_nIDs = len(editedIDs) + len(dat.wdID.unique()) + len(wo_lot.wdID.unique())
    editedIDs1 = read_text_file(file = outpath+f"\\matched\\{setID}_edited_1.txt")
    final_gdf.loc[final_gdf.wdID.isin(editedIDs1), 'code'] = 1
    if export:
        final_gdf.to_file(os.path.join(inpath, "output", "final", f"mapped_wd_{setID}.shp"), index=False)
    return final_gdf, toCheck, matched_gdf, digitized_nIDs, unmatchedIDs, issueIDs

def run_Tier3_4_final(setID, nm_to_add):
    """
    combine matched and digitized records and review
    gdf: the final shapefile that combined both automatic matches and digitized records
    """
    start = time.time()
    revpath = inpath + f'\GIS\ArcGIS Pro Project\DataReview\{setID}.gdb'
    wd = combine_wd_tables(setID=setID, nm_to_add=nm_to_add)
    matched = gpd.read_file(inpath + f'\\output\matched\matched_records_{setID}_edited.shp')
    mapped0 = [lyr for lyr in fiona.listlayers(revpath) if (lyr not in [f'{setID}_wo_lot', f'{setID}_partial']) and ('L' not in lyr)]
    partial = gpd.read_file(revpath, layer=f'{setID}_partial')
    edited = read_text_file(outpath+f'\\matched\\{setID}_edited.txt')
    editedIDs1 = read_text_file(file = outpath+f"\\matched\\{setID}_edited_1.txt")
    edited = edited + editedIDs1
    gdf, toCheck, matched_gdf, digitized_nIDs, unmatchedIDs, issueIDs = combine_matched_digitized(setID=setID, editedIDs=edited, nm_to_add=nm_to_add)
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
        gdf.to_file(inpath + '\\test\point.shp')
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
    res = str(cnt_dict[cnty]).zfill(2)
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
        lots_to_compare = ttdf.loc[ttdf.trsqq==trsqq_to_compare, 'ORTaxlot'].values
        lots_to_compare = list(map(get_lot_number_from_taxlot, lots_to_compare))
        if trsqq_to_compare_c == trsqq_to_check_c:
            print("trsqq matched, checking county code...")
            cnty_code = int(get_county_code_from_lonlat(lon, lat))
            county_to_compare = [key for key, value in cnt_dict.items() if value == cnty_code][0]
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
        r1_df.to_csv(os.path.join(inpath + '\\output\\to_review\\', f'unmatched_df_{setID}_r1_{ml}.csv'), index=False)
        r2_df.to_csv(os.path.join(inpath + '\\output\\to_review\\', f'unmatched_df_{setID}_r2_{ml}.csv'), index=False)
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
        taxlots_to_review_2.to_file(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_{ml}.shp'))
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
    df = df[df.county.isin(OR_counties)]
    outdf = df.copy()[['wetdet_delin_number', 'trsqq', 'parcel_id', 'county', 'latitude', 'longitude', 'DecisionLink', 'record_ID', 'IDyear']]
    outdf.loc[:,'correct_type'], outdf.loc[:,'correction'] = zip(*outdf.apply(lambda row: review_wd_record_w_coord(wd_id = row.wetdet_delin_number, county_to_check = row.county, trsqq_to_check = row.trsqq, parcel_IDs_to_check = row.parcel_id, lon = row.longitude,lat = row.latitude,taxlot = taxlot, year = row.IDyear), axis = 1))
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
        outdf.to_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_{ml}_0.csv'), index=False)
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
    notes = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'unmatched_df_{setID}_{s}_{ml}_notes.csv'))
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
        df.to_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_{s}_{ml}_1.csv'), index=False)
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
    rev_df = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_{ml}_0.csv'))
    nt_df = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'unmatched_df_{setID}_r2_{ml}_notes.csv'))
    df = df.copy()[~df.wetdet_delin_number.isin(nt_df.wetdet_delin_number.unique())]
    rev_df = rev_df.copy()[~rev_df.wetdet_delin_number.isin(nt_df.wetdet_delin_number.unique())]
    rev_df = rev_df[['wetdet_delin_number', 'cor_trsqq']]
    ndf = df.merge(rev_df, on='wetdet_delin_number')
    selectedID = ndf.cor_trsqq.astype(str) != 'nan'
    ndf.loc[selectedID, 'trsqq'] = ndf.loc[selectedID, 'cor_trsqq'].apply(lambda x: x.rstrip('0'))
    ndf.drop(columns='cor_trsqq', inplace=True)
    rev_df = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_{ml}_1.csv'))
    rdf = ndf.append(rev_df, ignore_index = True)
    if export:
        rdf.to_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_{ml}_2.csv'), index=False)
    return rdf

def combine_corrected_unmatched(setID, ml, skip=True, export=True):
    """
    combine the output from correct_unmatched and update_unmatched_df_r2
    setID: the setID of the unmatched records
    ml (missing lot): whether the unmatched records are missing parcel id in digit
    """
    rev_df1 = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r1_{ml}_1.csv'))
    # if skip update_unmatched_df_r2
    if skip:
        rev_df2 = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_{ml}_1.csv'))
    else:    
        rev_df2 = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_{ml}_2.csv'))
    df = rev_df1.append(rev_df2, ignore_index = True)
    if export:
        df.to_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_{ml}.csv'), index=False)
    return df

def review_WD_record_via_Pro(gdf, wdID):
    """
    review the unmatched records via Pro
    gdf: the output from split_unmatched_df
    wdID: the wetdet_delin_number of the unmatched records
    return the gdf of the unmatched records
    """
    gdf = gdf[gdf.wdID == wdID]
    gdf.to_file(os.path.join(inpath, 'output', 'wd_shp', wdID + '.shp'))
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
    taxlots_to_review_2.to_file(os.path.join(inpath, 'output', 'to_review', f'review_unmatched_{setID}_r2_{ml}.shp'))
    return taxlots_to_review_2

def adjust_taxlot(tx, ty):
    """
    adjust the taxlot with correct trsqq and taxlot with sheet number
    """
    res = ty
    if tx in tsq_dst:
        pattern = ty.split('--')[1][1:]
        tid_list = ttdf[ttdf.trsqq == tx].ORTaxlot.unique()
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
    if 'N' not in unmatched_df.missinglot.unique():
        print("All unmatched records are without taxlots, skip running Tier 2 Step 1...")
        return None, None
    else:
        r1_df, r2_df = split_unmatched_df(unmatched_df, ml='N', setID=setID)
        rev_r2 = review_unmatched_df_r2(r2_df, all_taxlot, setID, ml='N', export=True)
        taxlots_to_review = get_taxlot_to_check_r2(rev_r2, all_taxlot, setID, ml='N')
        return r1_df, r2_df

def run_Tier2_step3(r1_df, r2_df, setID, nm_to_add, wd, all_taxlot):
    """
    update the match
    """
    cor_r1 = correct_unmatched(r1_df, setID, s='r1', ml='N', export=True)
    cor_r2 = correct_unmatched(r2_df, setID, s='r2', ml='N', export=True)
    df = combine_corrected_unmatched(setID, ml='N')
    rev_df = reindex_data(df)
    matched = match_wd_data_with_taxlot(rev_df, setID, all_taxlot, export=True, update=True)
    unmatched_df = report_unmatched(matched, setID, nm_to_add, mute = False)
    matched_toReview = matched[matched.notes.notnull()] 
    wd_toReview = wd[wd.wetdet_delin_number.isin(matched_toReview.wdID.unique())]
    wd_toReview.to_csv(outpath + f'\\to_review\\partial_matched_{setID}.csv', index=False)
    unmatched_df.to_csv(os.path.join(inpath + '\\output\\to_review\\', f'unmatched_df_{setID}_2.csv'), index=False)
    return matched, unmatched_df

def report2DSL(setID):
    """
    generate the correction report for DSL
    """
    r2_0 = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_N_0.csv'))
    r1_notes = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'unmatched_df_{setID}_r1_N_notes.csv'))
    r2_notes = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'unmatched_df_{setID}_r2_N_notes.csv'))
    cordf = r1_notes.append(r2_notes, ignore_index=True)
    matched = gpd.read_file(outpath + f'\\matched\\matched_records_{setID}.shp')
    corrected = matched[matched.record_ID.isin(cordf.record_ID.values)][['wdID','trsqq', 'parcel_id', 'record_ID']].drop_duplicates(ignore_index=True)
    report = cordf.merge(corrected[['trsqq', 'parcel_id', 'record_ID']], on='record_ID')
    report['trsqq'] = report.trsqq.apply(lambda x: x.rstrip('0'))
    report.to_csv(outpath + f'\\corrected\\corrected_{setID}.csv', index=False)
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

def unique(list1):
    """
    This function takes a list and returns a list of unique values
    """
    x = np.array(list1)
    return list(np.unique(x))

def get_lot_numbers(x):
    """
    This function takes a string or integer and returns a list of lot numbers
    function to get all the lot numbers
    x is parcel_id
    """
    #print(x)
    if x is None:
        res = None
    elif type(x) is int:
        s = str(x)
        if len(str(x)) > 5:
            idx = [i for i, char in enumerate(s) if char != '0']
            lot_list = []
            for i in range(len(idx)-1):
                lot_list.append(s[idx[i]:idx[i+1]])
            lot_list.append(s[idx[len(idx)-1]:])
            res = lot_list
        else:
            res = [s]
    else:
        # remove parenthesis from text
        if '(' in str(x):
            txt = x.replace('(','').replace(')','')
        else:
            txt = str(x)
        # split the text
        lot_list = []
        for r in re.split(",|, | ", txt):
            if '-' in r:
                start, end = r.split('-')
                if start.isdigit() and end.isdigit():
                    lot_list += list(range(int(start), int(end)+1))
                    lot_list = list(map(lambda x: str(x), lot_list))
            else:
                lot_list.append(r)
        # remove text elements
        l = []
        # in case there are still number-letter strings (e.g., '1a')
        for t in [lot for lot in lot_list if ~lot.isnumeric()]:
            if any(c.isdigit() for c in t):
                l.append(re.sub('\D', '', t))
        res = unique([lot for lot in lot_list if lot.isnumeric()] + l)
        r = re.search('ROW|RR', x, re.IGNORECASE)
        w = re.search('Water', x, re.IGNORECASE)
        l = re.search('RAIL', x, re.IGNORECASE)
        for srch, msg in zip([r, w, l], ['ROADS','WATER','RAILS']):
            if srch:
                #res += msg
                res.append(msg)
    return res

def check_duplicates(v):
    """
    check duplicates in a list
    """
    itemlist = [item for item, count in collections.Counter(v).items() if count > 1]
    return itemlist, len(itemlist)

def read_wd_table(setID, file):
    """
    read wd tables, recordID is used for single tables
    """
    datafile = os.path.join(wdpath, setID, file)
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
    
def get_tr_code(x, code='t'):
    """
    get township or range code from the taxlot
    """
    nms = re.findall('\d+', x)
    if code=='t':
        nm1 = nms[0][0:2]
    else:
        nm1 = nms[1][0:2]
    lts = re.findall("[a-zA-Z]+", x)
    if code=='t':
        lts2 = lts[0]
    else:
        dirlst = ["E", "W", "S", "N"]
        if lts[1] in dirlst:
            lts2 = lts[1]
        else:
            lts2 = sample(dirlst, 1)[0]
    
    if len(nm1) == 1:
        tr1 = '0' + nm1
    elif len(nm1) == 3:
        tr1 = nm1[1:3]
    else:
        tr1 = nm1
    
    if ('V' in lts2) or ('Y' in lts2):
        tr2 = '.50'
        tr3 = lts2[1]
    elif('X' in lts2) or ('Z' in lts2):
        if any([x in lts2 for x in ['XS', 'ZN', 'XE', 'ZW']]):
            tr2 = '.75'
        else:
            tr2 = '.25'
    else:
        tr2 = '.00'
        tr3 = lts2

    res = tr1 + tr2 + tr3
    return res
   
def get_s_code(x):
    """
    get section code from trsqq code
    """
    if len(x) <= 7:
        s = '00'
    else:
        nms = re.findall('\d+', x)
        n = len(nms)
        k = len(nms[1])
        if(n < 3) & (k > 2):
            nm1 = nms[1][(k-2):(k+1)]
        else:
            nm1 = nms[2]
        if len(nm1) == 1:
            s = '0' + nm1
        else:
            s = nm1
    return s       

def get_qq_code(x):
    """
    get QQ code from trsqq code
    """
    dirlst = ["E", "W", "S", "N"]
    nms = re.findall('\d+', x)
    lts = re.findall("[a-zA-Z]+", x)
    if (len(lts) == 2) & (lts[1] not in dirlst):
        if len(lts[1]) == 2:
            qq = lts[1]
        else:
            qq = lts[1] + '0'
    elif len(nms[2]) > 2:
        if nms[2] == 4:
            qq = nms[2][2:4]
        else:
            qq = nms[2][2] + '0'
    elif len(lts) == 3:
        if len(lts[2]) == 2:
            qq = lts[2]
        else:
            qq = lts[2] + '0'
    else:
        if len(nms[0]) == 1:
            t = '0' + nms[0]
        else:
            t = nms[0]
        if len(nms[1]) == 1:
            r = '0' + nms[1]
        else:
            r = nms[1]
        if len(nms[2]) == 1:
            s = '0' + nms[2]
        else:
            s = nms[2]
        trsqq = t + lts[0] + r + lts[1] + s
        qq =  '{:0<10}'.format(trsqq)[8] + '{:0<10}'.format(trsqq)[9]
    return qq

def convert_trsqq(x):
    """
    convert township, range, section, and quarter-quarter to the taxlot id format
    """
    #print(x)
    x = '{:<08s}'.format(x)
    #x = re.sub("V|X|Y|Z", "", x)
    xt = get_tr_code(x, code='t') + get_tr_code(x, code='r') + get_s_code(x) + get_qq_code(x)
    return xt[:16]

def create_ORTaxlot(cnt_code, trsqq, lot):
    """
    create the taxlot id based on the county code, township, range, section, and lot number
    """
    #print(f'County {cnt_code}, TRSQQ {trsqq}, Lot {lot}')
    part1 = str(int(cnt_code)).zfill(2) + convert_trsqq(trsqq) 
    taxlotID = part1 + '--' + ('000000000' + lot)[-9:]
    tid_dst_2 = [x[-len(lot):] for x in tid_dst_1]    
    if (taxlotID not in all_txid):
        if (part1 in tid_dst_0) and (lot in tid_dst_2):
            txl_ID_lst = [tid for tid in tid_dst if (re.search(part1, tid)) and (re.search(lot, tid))]
            if len(txl_ID_lst) > 0:
                taxlotID = txl_ID_lst[0]
            else:
                ntaxlotID = taxlotID[:6]+taxlotID[7:12]+taxlotID[13:18]+'00--'+taxlotID.split('--')[1]
                if ntaxlotID in all_txid:
                    taxlotID = ntaxlotID         
    return taxlotID

def reindex_data(wd_dt):
    """
    reindex the data based on the number of lots in each parcel id
    wd_dt is from read_wd_table or reorganize_tocheck
    make sure the records are with lot numbers
    return reindexed data
    """
    # get a list of lot numbers in each parcel id record
    # will need to review the records without any parcel ids
    selectedID = wd_dt.parcel_id.astype(str) != 'nan'
    wd_dt = wd_dt.copy()[selectedID]
    wd_dt.loc[:, 'lots'] = wd_dt['parcel_id'].apply(lambda x: get_lot_numbers(x))
    # repeat the rows based on the number of lot numbers
    ndf = wd_dt.reindex(wd_dt.index.repeat(wd_dt.lots.str.len()))
    # add the column to list the lot for all
    ndf.loc[:, 'lot'] = list(chain.from_iterable(wd_dt.lots.values.tolist()))
    ndf.loc[:, 'lots'] = ndf['lots'].apply(lambda x: ', '.join(dict.fromkeys(x).keys()))
    # get county code
    ndf.loc[:, 'cnt_code'] = ndf.county.map(cnt_dict)
    # get OR taxlot IDs for wd data
    ndf = ndf[~ndf.cnt_code.isnull()]
    ndf.loc[:, 'ORTaxlot'] = ndf[['cnt_code', 'trsqq', 'lot']].apply(lambda row: create_ORTaxlot(cnt_code=row.cnt_code, trsqq=row.trsqq, lot=row.lot), axis = 1)
    return ndf

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


def clean_wd_table(setID, file):
    """
    function to clean up wd data by single file
    """
    start = time.time()
    wd_dt = read_wd_table(setID, file)
    # this will help identify problematic records with numbers
    selectedID = wd_dt.parcel_id.astype(str) != 'nan'
    wd_dt.loc[selectedID, 'notes'] = wd_dt.copy()[selectedID]['parcel_id'].apply(lambda x: make_notes(x))
    wd_dt['county'] = wd_dt['county'].apply(lambda x: x.title())
    wd_dt = wd_dt[wd_dt.county.isin(OR_counties)]
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
        ndf = ndf[ndf.county.isin(OR_counties)]
        end = time.time()
        res = ndf
        #print(f'cleaned up wd data in {file} and it took about {end - start} seconds')
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

def combine_wd_tables(setID, nm_to_add, raw=True):
    """
    Combine all the wd tables in the set to review unique records, record_ID is used for combined tables
    use this function when reindex is not neccessary
    nm_to_add is the number of previous records (records from the previous sets; same to all functions with this variable)
    """  
    if raw:
        frames = []
        files = list_files(os.path.join(wdpath, setID))
        # in case there are unidentified files
        files = [file for file in files if '~$' not in file]
        for file in files:
            datafile = os.path.join(wdpath, setID, file)
            xl = pd.ExcelFile(datafile)
            wd_dt = pd.read_excel(datafile, sheet_name=xl.sheet_names[1])
            frames.append(wd_dt)
        wd_df = pd.concat(frames, ignore_index=True)
        # this creates unique IDs for all the records in the same set
        wd_df.loc[:, 'record_ID'] = range(1, wd_df.shape[0] + 1) 
        wd_df.loc[:, 'record_ID'] = wd_df.copy().loc[:, 'record_ID'] + nm_to_add
        selectedID = wd_df.parcel_id.astype(str) != 'nan'
        wd_df.loc[selectedID, 'notes'] = wd_df[selectedID]['parcel_id'].apply(lambda x: make_notes(x))
        wd_df.loc[:, 'recyear'] = wd_df.received_date.apply(lambda x: x.year)
        # get year from the wd ID 
        wd_df.loc[:, 'IDyear'] = wd_df.wetdet_delin_number.apply(lambda x: x[2:6]) 
        selectedID = wd_df.parcel_id.astype(str) != 'nan'
        wd_df.loc[selectedID, 'missinglot'] = wd_df[selectedID].parcel_id.apply(lambda x: without_lots(x))
    else:
        wd_df = pd.read_csv(wdpath+f'\\Corrected_by_Set\\{setID}.csv')
    
    return wd_df

def read_taxlot(year, mute=True):
    """
    function to read taxlot data
    """
    start = time.time()
    txfilepath = os.path.join(txpath, 'Taxlots' + str(year) + '.gdb')
    if year < 2016:
        tx_dt = gpd.read_file(txfilepath, layer=f'Taxlots{year}')
    else: 
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
                   yearstart=2016, 
                   yearend=2023,
                   skips=[2010, 2013]):
    """
    combine taxlots from all years
    """
    frames = []
    for year in range(yearstart, yearend):
        if year not in skips:
            print(year)
            tx_dt = read_taxlot(year)
            if 'Year' not in tx_dt.columns:
                tx_dt['year'] = str(year)
            else:
                tx_dt.rename(columns={'Year': 'year'}, inplace=True)
            frames.append(tx_dt[['year', 'ORTaxlot', 'geometry']])
    df = pd.concat(frames, ignore_index=True)
    gdf = gpd.GeoDataFrame(df, crs="EPSG:2992", geometry='geometry')
    if exportID:
        with open(os.path.join(inpath, "ORTaxlot.pkl"), "wb") as f:
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
    files = list_files(os.path.join(wdpath, setID))
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
    found = [txid for txid in tocheck_txid if txid in all_txid]
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
            matched = gpd.read_file(os.path.join(inpath + f'\\{outfolder}\\', f'matched_records_{setID}.shp'))
            ngdf = matched.append(ngdf[selcols], ignore_index = True)
            ngdf = ngdf[~ngdf.geometry.isnull()]
        if export: 
            # the update will overwrite the first output
            ngdf[selcols].to_file(os.path.join(inpath + f'\\{outfolder}\\', f'matched_records_{setID}.shp'), driver='ESRI Shapefile')  
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
    matched_rID = gdf[~gdf.geometry.isnull()].record_ID.unique()
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
        sorted_df.to_csv(os.path.join(inpath + '\\output\\to_review\\', f'unmatched_df_{setID}.csv'), index=False)  
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
    setgdf = gpd.read_file(os.path.join(inpath, 'GIS', 'Join_Statewide.gdb'), layer=f'WD_{setID}_Combined')
    setgdf.loc[setgdf.Record_ID.astype(str) != 'nan', 'Record_ID'] = setgdf[setgdf.Record_ID.astype(str) != 'nan'].Record_ID.astype('int64', copy=False)
    matched_rID = gdf.record_ID.unique()
    # missed IDs in the existing data that is not nan
    missed_gdf = setgdf[setgdf.Record_ID.astype(str) != 'nan'][~setgdf.Record_ID.isin(matched_rID)]
    missedID = missed_gdf.Record_ID.unique()
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
            missed_gdf[~missed_gdf.geometry.isnull()][selcols].to_file(os.path.join(inpath + f'\\{outfolder}\\',
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
            added_gdf.to_file(os.path.join(inpath + f'\\{outfolder}\\', f'added_records_in_{setID}_res.shp'), 
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
    setgdf = gpd.read_file(os.path.join(inpath, 'GIS', 'Join_Statewide.gdb'), layer=f'WD_{setID}_Combined')
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
        mapped[selcols].to_file(os.path.join(inpath + f'\\{outfolder}\\', f'mapped_in_{setID}.shp'), 
                                                  driver='ESRI Shapefile')
    return mapped

def check_corrected_data(df, setID, all_taxlot, nm_to_add, export=False):
    """
    update match on the corrected data
    input df is reindexed from combined_reindexed_data
    return matched data with corrected info
    """ 
    corrected = pd.read_excel(os.path.join(inpath, 'DSL data originals', 'Data corrections feedback to DSL', 
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
        ngdf[~ngdf.geometry.isnull()][selcols].to_file(os.path.join(inpath + '\\output\\matched\\', f'matched_records_{setID}.shp'), driver='ESRI Shapefile')
    return ndf, cor_df, comIDs, df_wlots_to_check

Tdict = dict(zip(['.25S', '.50S', '.50N', '.75S'], ['Z', 'V', 'Y', 'X']))
Rdict = dict(zip(['.25E', '.50E', '.50W', '.75E'], ['Z', 'V', 'Y', 'X']))

def taxlot2trsqq(x):
    """
    convert taxlot to trsqq
    """ 
    #print(x)
    rcs = ['.25', '.50', '.75'] # rcs = ratio codes
    ptns = '.25|.50|.75'
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
    taxlotIDs_to_search = all_txid
    taxlotIDs_cleaned = unique([txID for txID in taxlotIDs_to_search if len(txID) == 29])
    trsqq = list(map(lambda x: taxlot2trsqq(x), taxlotIDs_cleaned))
    trsqq_dict = dict(zip(trsqq, taxlotIDs_cleaned))
    df = pd.DataFrame({'ORTaxlot':taxlotIDs_cleaned, 'trsqq':trsqq})
    with open(os.path.join(inpath, "trsqq_list.pkl"), "wb") as f:
        pickle.dump(trsqq, f)
    with open(os.path.join(inpath, "trsqq_dict.pkl"), "wb") as f:
        pickle.dump(trsqq_dict, f)
    df.to_csv(os.path.join(inpath, "trsqq_df.csv"), index=False)
    return trsqq, trsqq_dict, df

def get_maybe_taxlot(trsqq_to_check):
    """
    get the maybe taxlot from the trsqq_to_check
    this is a test function
    """    
    closematch = difflib.get_close_matches(trsqq_to_check, trsqq)
    trsqq_matched = unique(closematch)
    checktaxlot = [*map(trsqq_dict.get, trsqq_matched)]
    string_to_search = trsqq_matched[0][0:8]
    search_res = unique([i for i in unique(ttdf.trsqq) if string_to_search in i])
    if len(checktaxlot) == 1 and len(search_res)==1:
        values = [i for i in trsqq_dict if trsqq_dict[i]==checktaxlot[0]]
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