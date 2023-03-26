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
from PyPDF2 import PdfFileReader, PdfFileWriter
import fiona
import webbrowser
import time

inpath = r'L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping'
wdpath = inpath + '\\DSL data originals'
txpath = inpath + '\\GIS\ORMAP_data\ORMAP_Taxlot_Years'
yearstart = 2016
yearend = 2023
#outfolder = 'test'
outfolder = 'output\\matched'
# create a spreadsheet to create a dictionary for the match between county name and code
cnt_ID = pd.read_excel(r'T:\DCProjects\EPA-WD\CNT_Code.xlsx')
# create a dictionary to look up county code
cnt_dict = dict(zip(cnt_ID.COUNTY, cnt_ID.ID))

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
tsq_dst = ttdf[ttdf.ORTaxlot.isin(tid_dst)].trsqq.unique()

revpath = inpath + '\GIS\ArcGIS Pro Project\DataReview\DataReview.gdb'
pdf_outpath = r'L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping\output\pdf'
with open(os.path.join(inpath, "ORTaxlot.pkl"), "rb") as f:
        all_txid = pickle.load(f)

pd.options.mode.chained_assignment = None

################################################ Tier 3 & 4 #####################################################
def review_loop_r1(df):
    for wdID in df.wetdet_delin_number.unique():
        print(wdID)
        print(check_unmatched_r1(wdID = wdID, df = df))
        user_input = input("Press 'p' to pause or any other key to continue...")
        if user_input == 'p':
            while True:
                user_input = input("Press 'c' to continue...")
                if user_input == 'c':
                    break
        time.sleep(1) # wait for 1 second between iterations

def check_unmatched_r1(wdID, df):
    url = df.loc[df.wetdet_delin_number == wdID, 'DecisionLink'].values[0]
    selcols = ['county', 'trsqq', 'parcel_id', 'latitude', 'longitude', 'record_ID', 'notes']
    if str(url) == 'nan':
        print('Decision link is not available')
    else:
        webbrowser.open(url)
    return df.loc[df.wetdet_delin_number == wdID, selcols]

def review_loop(df):
    df = df.reset_index()
    for i in range(df.shape[0]):
        wdID = df.loc[i, 'wetdet_delin_number']
        print(wdID)
        print(check_review_notes_r2n(wdID = wdID, df = df))
        user_input = input("Press 'p' to pause or any other key to continue...")
        if user_input == 'p':
            while True:
                user_input = input("Press 'c' to continue...")
                if user_input == 'c':
                    break
        time.sleep(1) # wait for 1 second between iterations

def check_review_notes_r2n(wdID, df):
    url = df.loc[df.wetdet_delin_number == wdID, 'DecisionLink'].values[0]
    if str(url) == 'nan':
        print('Decision link is not available')
    else:
        webbrowser.open(url)
    return df.loc[df.wetdet_delin_number == wdID, ['correct_type', 'correction', 'cor_trsqq']]

def check_completeness(setID='003', a=3):
    partial = gpd.read_file(revpath, layer=f'Set{setID}_partial')
    mapped1 = list(partial.wdID.unique())
    mapped0 = [lyr for lyr in fiona.listlayers(revpath) if (lyr not in [f'Set{setID}_wo_lot', f'Set{setID}_partial']) and ('L' not in lyr)]
    mapped2 = list(map(lambda x: x.replace('_', '-'), mapped0))
    mapped = mapped1 + mapped2
    matched = gpd.read_file(inpath + f'\\output\matched\matched_records_Set{setID}_edited.shp')
    pct = (len(mapped)+a)/len(matched[~matched.notes.isnull()].wdID.unique())
    print(f'{round(pct*100, 1)}% completed...')
    return sorted(mapped), len(mapped)

def extract_page_from_locPath(filePath, pageNm, wdID):
    pdf_file = PdfFileReader(filePath)
    pageObj = pdf_file.getPage(pageNm)
    pdf_writer = PdfFileWriter()
    pdf_writer.addPage(pageObj)
    output = f'{pdf_outpath}\\{wdID}_{pageNm+1}.pdf'
    with open(output, 'wb') as output_pdf:
        pdf_writer.write(output_pdf) 

def extract_page_from_docLink(url, pageNm, wdID):
    response = requests.get(url=url, timeout=120)
    on_fly_mem_obj = io.BytesIO(response.content)
    pdf_file = PdfFileReader(on_fly_mem_obj)
    pageObj = pdf_file.getPage(pageNm)
    pdf_writer = PdfFileWriter()
    pdf_writer.addPage(pageObj)
    output = f'{pdf_outpath}\\{wdID}_{pageNm+1}.pdf'
    with open(output, 'wb') as output_pdf:
        pdf_writer.write(output_pdf) 

def review_mapped(setID):
    revpath = inpath + f'\GIS\ArcGIS Pro Project\DataReview\{setID}.gdb'
    mapped0 = [lyr for lyr in fiona.listlayers(revpath) if (lyr not in [f'{setID}_wo_lot', f'{setID}_partial']) and ('L' not in lyr)]
    for wID in mapped0:
        gdf = gpd.read_file(revpath, layer=wID)
        if 'wdID' in gdf.columns:
            if len(gdf.wdID.unique()) > 1:
                print(wID)

def revise_single_partial_file(wID):
    gdf = gpd.read_file(revpath, layer=wID)
    gdf = gdf.to_crs(epsg=2992)
    selcols = ['Shape_Length', 'Shape_Area']
    if 'wdID' in gdf.columns:
        df = gdf.copy()[selcols +  ['wdID']]
    else:
        df = gdf.copy()[selcols]
        df.loc[:,'wdID'] = wID.replace('_', '-')
    df.loc[:,'geometry'] = gdf.loc[:,'geometry']
    return df

def merge_single_partial_file(wIDlist):
    df = pd.DataFrame()
    for wID in wIDlist:
        df=pd.concat([df, revise_single_partial_file(wID)], ignore_index=True)
    gdf = gpd.GeoDataFrame(df, crs="EPSG:2992", geometry='geometry')
    return gdf

# require the edited matched records, digitized partial taxlots, taxlots without lot IDs, and the list of issue IDs
def combine_matched_digitized(setID, editedIDs, nm_to_add, export=True):
    mapped0 = [lyr for lyr in fiona.listlayers(revpath) if (lyr not in [f'{setID}_wo_lot', f'{setID}_partial']) and ('L' not in lyr)]
    matched = gpd.read_file(inpath + f'\\output\matched\matched_records_{setID}_edited.shp')
    partial = gpd.read_file(revpath, layer=f'{setID}_partial')
    partial = partial.to_crs(epsg=2992)
    gdf = merge_single_partial_file(mapped0)
    dat = gdf.append(partial, ignore_index=True)
    edited_gdf = matched[matched.wdID.isin(editedIDs)]
    edited_gdf = edited_gdf[['wdID', 'geometry']].dissolve('wdID')
    edited_gdf.loc[:, 'wdID'] = edited_gdf.index
    data1 = edited_gdf.append(dat[['wdID', 'geometry']], ignore_index=True)
    matched_df = matched[['wdID', 'geometry']].dissolve('wdID')
    matched_df.loc[:, 'wdID'] = matched_df.index
    excluded = [wdID for wdID in data1.wdID.unique() if wdID in matched_df.wdID.unique()]
    issues = pd.read_csv(os.path.join(inpath, "output", "to_review", f"{setID}_Mapping_Issues.csv"))
    issueIDs = list(issues.wetdet_delin_number.unique())
    matched_gdf = matched_df[~matched_df.wdID.isin(excluded+issueIDs)]
    data2 = matched_gdf.append(data1, ignore_index=True)
    wo_lot = gpd.read_file(revpath, layer=f'{setID}_wo_lot')
    wo_lot = wo_lot.to_crs(epsg=2992)
    wd = combine_wd_tables(setID=setID, nm_to_add=nm_to_add)
    data3 = data2.append(wo_lot[['wdID', 'geometry']], ignore_index=True)
    unmatchedIDs = [wdID for wdID in wd.wetdet_delin_number.unique() if wdID not in data3.wdID.unique()]
    toCheck = [ID for ID in unmatchedIDs if ID not in issueIDs]
    digitized_nIDs = len(editedIDs) + len(dat.wdID.unique()) + len(wo_lot.wdID.unique())
    if export:
        data3.to_file(os.path.join(inpath, "output", "final", f"mapped_wd_{setID}.shp"))
    return data3, toCheck, digitized_nIDs, issueIDs

################################################ Tier 2 #####################################################
def get_point_from_lonlat(lon, lat, export=True):
    df = pd.DataFrame([[lon, lat]], columns=['Longitude', 'Latitude'])
    gdf = gpd.GeoDataFrame(df, crs="EPSG:4326", geometry=gpd.points_from_xy(df.Longitude, df.Latitude))
    gdf = gdf.to_crs(epsg=2992)
    if export:
        gdf.to_file(inpath + '\\test\point.shp')
    return gdf

# point in polygon - WD point in taxtlot
# require taxtlot
def extract_taxlot_info(wd_pt, taxlot, year):
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

def get_county_code_from_lonlat(lon, lat):
    geolocator = Nominatim(user_agent="geoapiExercises")
    location = geolocator.reverse(str(lat)+","+str(lon))
    address = location.raw['address']
    county = address.get('county', '')
    cnty = county.split(' ')[0]
    res = str(cnt_dict[cnty]).zfill(2)
    return res

def separate_numbers_letters(x):
    numbers = re.findall('\d+', x)
    letters = re.findall("[a-zA-Z]+", x)
    return numbers, letters

def find_different_indices(list1, list2):
    different_indices = []
    for i, (a, b) in enumerate(zip(list1, list2)):
        if a != b:
            different_indices.append(i)
    return different_indices

def pad_string(string, length=10):
    if len(string) < length:
        padded_string = string.ljust(length, '0')
    else:
        padded_string = string
    return padded_string

def combine_lists(list1, list2):
    combined_list = list(zip(list1, list2))
    return combined_list

def remove_tuple_format(input_list):
    output_list = []
    for item in input_list:
        if isinstance(item, tuple):
            output_list.extend(item)
        else:
            output_list.append(item)
    return output_list

def get_list_elements_by_index(input_list, index_list):
    output_list = [input_list[i] for i in index_list]
    return output_list

def split_trsqq(trsqq_to_check):
    numbers1, letters1 = separate_numbers_letters(trsqq_to_check[:-2])
    letters1.append(trsqq_to_check[-2:])
    trsqq_to_check_lst = remove_tuple_format(combine_lists(numbers1, letters1))
    return trsqq_to_check_lst

def compare_trsqq(trsqq_to_check, trsqq_to_compare):
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
    delimiter1 = ', '
    delimiter2 = ' and '
    position1 = 1
    position2 = -1
    res = delimiter1.join(my_list[:position1]) + delimiter1 + delimiter1.join(my_list[position1:position2]) + delimiter2 + my_list[position2]
    return res

trsqq_correction_dict = dict(zip(list(range(0, 6)), ['township number', 'township direction', 'range number', 'range direction', 'section number', 'QQ']))

def report_trsqq_correction(trsqq_to_check, trsqq_to_compare, to_correct=False):
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
    #print(trsqq_to_check)
    wd_pt  = get_point_from_lonlat(lon = lon, lat = lat)
    tID = extract_taxlot_info(wd_pt = wd_pt, taxlot = taxlot, year = year)
    trsqq_to_check_c = pad_string(trsqq_to_check)
    trsqq_to_compare = trsqq_dict[tID]
    trsqq_to_compare_c = pad_string(trsqq_to_compare)
    res = report_trsqq_correction(trsqq_to_check_c, trsqq_to_compare_c, to_correct=True)    
    return res, tID

# this function only works when the coordinates are accurate and one-on-one match among WD ID, trsqq, parcel IDs, and coordinate
# limitation - one WD record (possible with multiple records with different trasqq and parcel IDs) has only one coordindate
def review_wd_record_w_coord(wd_id, county_to_check, trsqq_to_check, parcel_IDs_to_check, lon, lat, taxlot, year):
        print(f'reviewing {wd_id}')
        wd_pt  = get_point_from_lonlat(lon = lon, lat = lat)
        tID = extract_taxlot_info(wd_pt = wd_pt, taxlot = taxlot, year = year)
        if "away" not in tID:
            trsqq_to_check_c = pad_string(trsqq_to_check)
            lots_to_check = get_lot_numbers(parcel_IDs_to_check)
            trsqq_to_compare = trsqq_dict[tID]
            trsqq_to_compare_c = pad_string(trsqq_to_compare)
            lots_to_compare = ttdf.loc[ttdf.trsqq==trsqq_to_compare, 'ORTaxlot'].values
            lots_to_compare = list(map(get_lot_number_from_taxlot, lots_to_compare))  
            if trsqq_to_compare_c == trsqq_to_check_c:
                print("trsqq matched, checking county code...")
                cnty_code = int(get_county_code_from_lonlat(lon, lat))
                county_to_compare = [key for key, value in cnt_dict.items() if value == cnty_code]
                # need to check the typos in the county name first
                if county_to_check == county_to_compare[0]:
                    print("county code is corrected, need to check lot numbers...")
                    if any([x not in lots_to_compare for x in lots_to_check]):
                        lots_to_correct = [x for x in lots_to_check if x not in lots_to_compare]
                        cor_type, cor_notes = "lot number", f'lot number {lots_to_correct} might be incorrect, the matched taxlot is {tID} for {trsqq_to_compare}'
                        print("lot numbers might be wrong...")
                    else:
                        notes = 'lot numbers seem to be correct, need to review'
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

# df is the output from report_unmatched
def split_unmatched_df(df, ml, setID, export=True):
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

# df is the second output from split_unmatched_df
def generate_taxlot_output(df, taxlot, setID, ml, export=False):
    df['pairs'] = list(zip(df['IDyear'].astype(str), df['ORTaxlot']))
    taxlots_to_review = taxlot[taxlot[['year', 'ORTaxlot']].apply(tuple, axis=1).isin(df.pairs.values)]
    taxlots_to_review_2 = taxlots_to_review.merge(df, on='ORTaxlot')
    taxlots_to_review_2.drop(columns=['pairs'], inplace=True)
    if export:
        taxlots_to_review_2.to_file(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_{ml}.shp'))
    return taxlots_to_review_2

def taxlot_from_coord(x):
    txid = x.split('the matched taxlot is ')[1].split(' for ')[0]
    trsqq = x.split(' for ')[1]
    return trsqq, txid

def trsqq_from_nearby_taxlot(x):
    txid = x.split(', about')[0].replace('coordinate might be incorrect, nearby taxlot is ', '')
    return trsqq_dict[txid], txid

# df is the second output from split_unmatched_df
def review_unmatched_df_r2(df, taxlot, setID, ml, export=True):
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
        outdf.to_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_{ml}_0.csv'), index=False)
    return outdf

trsqq_cor_dict = {'township number': 0, 
                  'township direction':2, 
                  'range number':3, 
                  'range direction':5,
                  'section number':6, 
                  'QQ':8}

def replace_str_index(text,index=0,replacement=''):
    return text[:index] + replacement + text[index+len(replacement):]

# df is the output from split_unmatched_df
# need to run review_unmatched_df_r2 and do some manual review work first to get the notes
# the document to review is review_unmatched_Set00?_r2_N_0.csv
def correct_unmatched(df, setID, s, ml, export=True):
    notes = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'unmatched_df_{setID}_{s}_{ml}_notes.csv'))
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
                    k = notes[sel2].shape[0]
                    if k > 1:
                        df.loc[sel, field] = df.loc[sel, field].apply(lambda x: pad_string(x))
                        cor_types = notes.loc[sel2, 'cor_type'].values
                        for cor_type in cor_types:
                            sel3 = sel2 & (notes.cor_type==cor_type)
                            ind = trsqq_cor_dict[cor_type]
                            df.loc[sel, field] = df.loc[sel, field].apply(lambda x: replace_str_index(x,
                                                                                                index=ind,
                                                                                                replacement=notes.loc[sel3, 'to'].values[0]))
                    else:
                        to = notes.loc[sel2, 'to'].values[0]
                        if to.isdigit():
                            to = to.zfill(2)
                        df.loc[sel, field] = df.loc[sel, field].apply(lambda x: x.replace(notes.loc[sel2, 'from'].values[0],
                                                                                    to))
        else:    
            sel = df.wetdet_delin_number == wdID
            fields = notes.loc[notes.wetdet_delin_number == wdID, 'field'].unique()
            for field in fields:
                sel2 = (notes.wetdet_delin_number==wdID) & (notes.field==field)
                k = notes[sel2].shape[0]
                if k > 1:
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
                    if to.isdigit():
                        to = to.zfill(2)
                    df.loc[sel, field] = df.loc[sel, field].apply(lambda x: x.replace(notes.loc[sel2, 'from'].values[0],
                                                                                to))
    if export:
        df.to_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_{s}_{ml}_1.csv'), index=False)
    return df

# need to run review_unmatched_df_r2 and do some manual work frist
# rev_df is the output from review_unmatched_df_r2, nt_df is from manual input
# df is the second output from split_unmatched_df
# this is a parallel run with correct_unmatched for r2_df
def update_unmatched_df_r2(df, setID, ml, export=True):
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

def combine_corrected_unmatched(setID,ml,export=True):
    rev_df1 = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r1_{ml}_1.csv'))
    rev_df2 = pd.read_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_r2_{ml}_2.csv'))
    df = rev_df1.append(rev_df2, ignore_index = True)
    if export:
        df.to_csv(os.path.join(inpath + '\\output\\to_review\\', f'review_unmatched_{setID}_{ml}.csv'), index=False)
    return df

def review_WD_record_via_Pro(gdf, wdID):
    gdf = gdf[gdf.wdID == wdID]
    gdf.to_file(os.path.join(inpath, 'output', 'wd_shp', wdID + '.shp'))
    return gdf

# df is the output from review_unmatched_df_r2
# make sure ORTaxlot is in the right format from df
def get_taxlot_to_check_r2(revdf, taxlot, setID, ml):
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
    df.loc[:, 'ORTaxlot'] = df.copy()[['trsqq', 'ORTaxlot']].apply(lambda row: adjust_taxlot(row.trsqq, row.ORTaxlot), axis=1)
    return df

################################################ Tier 1 #####################################################
# gdf below generally refers to the matched records

# function to list files or folders
def list_files(path, folder=False):
    f = []
    for (dirpath, dirnames, filenames) in walk(path):
        f.extend(filenames)
        break
    if folder:
        f = [x[0] for x in os.walk(path)]
    return f

# function to remove duplicates in the parcel ids
def unique(list1):
    x = np.array(list1)
    return list(np.unique(x))

# function to get all the lot numbers
# x is parcel_id
def get_lot_numbers(x):
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

# function to review duplicated records
# return the list and the list length
def check_duplicates(v):
    itemlist = [item for item, count in collections.Counter(v).items() if count > 1]
    return itemlist, len(itemlist)

# read wd tables, recordID is used for single tables
def read_wd_table(setID, file):
    datafile = os.path.join(wdpath, setID, file)
    xl = pd.ExcelFile(datafile)
    wd_dt = pd.read_excel(datafile, sheet_name=xl.sheet_names[1])
    wd_dt.loc[:, 'county'] = wd_dt.county.apply(lambda x: x.capitalize())
    # this will show the records in the original single table and the ID will be updated when the tables are combined
    wd_dt.loc[:, 'recordID'] = range(1, wd_dt.shape[0] + 1)
    return wd_dt

# get township or range code
def get_tr_code(x, code='t'):
    nms = re.findall('\d+', x)
    if code=='t':
        nm1 = nms[0]
    else:
        nm1 = nms[1]
    lts = re.findall("[a-zA-Z]+", x)
    if code=='t':
        lts2 = lts[0]
    else:
        lts2 = lts[1]
    
    if len(nm1) == 1:
        tr1 = '0' + nm1
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
   
# get section code
def get_s_code(x):
    if len(x) <= 7:
        s = '00'
    else:
        nms = re.findall('\d+', x)
        nm1 = nms[2]
        if len(nm1) == 1:
            s = '0' + nm1
        else:
            s = nm1
    return s       

def convert_trsqq(x):
    s = re.sub("V|X|Y|Z", "", x)
    xt = get_tr_code(x, code='t') + get_tr_code(x, code='r') + get_s_code(x) + '{:0<10}'.format(s)[8] + '{:0<10}'.format(s)[9]
    return xt[:16]

def create_ORTaxlot(cnt_code, trsqq, lot):
    taxlotID = str(cnt_code).zfill(2) + convert_trsqq(trsqq) + '--' + ('000000000' + lot)[-9:]
    return taxlotID

# wd_dt is from read_wd_table or reorganize_tocheck
# return reindexed data
def reindex_data(wd_dt):
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
    ndf.loc[:, 'ORTaxlot'] = ndf[['cnt_code', 'trsqq', 'lot']].apply(lambda row: create_ORTaxlot(cnt_code=row.cnt_code, trsqq=row.trsqq, lot=row.lot), axis = 1)
    return ndf

# add notes to the wd records
def make_notes(text):
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

# function to clean up wd data by single file
def clean_wd_table(setID, file):
    start = time.time()
    wd_dt = read_wd_table(setID, file)
    # this will help identify problematic records with numbers
    selectedID = wd_dt.parcel_id.astype(str) != 'nan'
    wd_dt.loc[selectedID, 'notes'] = wd_dt.copy()[selectedID]['parcel_id'].apply(lambda x: make_notes(x))
    ndf = reindex_data(wd_dt)
    # get year from the receive date
    ndf.loc[:, 'recyear'] = ndf.copy().received_date.apply(lambda x: x.year)
    # get year from the wd ID 
    ndf.loc[:, 'IDyear'] = ndf.wetdet_delin_number.apply(lambda x: x[2:6])
    ndf.loc[selectedID, 'missinglot'] = ndf.loc[selectedID, 'parcel_id'].apply(lambda x: without_lots(x))
    ndf['response_date'] = ndf['response_date'].dt.strftime("%Y-%m-%d")
    ndf['received_date'] = ndf['received_date'].dt.strftime("%Y-%m-%d")
    end = time.time()
    #print(f'cleaned up wd data in {file} and it took about {end - start} seconds')
    return ndf

# check whether the parcel ID is without lots
def without_lots(text):
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

# combine all the wd tables in the set to review unique records, record_ID is used for combined tables
# use this function when reindex is not neccessary
# nm_to_add is the number of previous records (records from the previous sets; same to all functions with this variable)
def combine_wd_tables(setID, nm_to_add):
    files = list_files(os.path.join(wdpath, setID))
    # in case there are unidentified files
    files = [file for file in files if '~$' not in file]
    frames = []
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
    # get year from the receive date
    wd_df.loc[:, 'recyear'] = wd_df.received_date.apply(lambda x: x.year)
    # get year from the wd ID 
    wd_df.loc[:, 'IDyear'] = wd_df.wetdet_delin_number.apply(lambda x: x[2:6]) 
    selectedID = wd_df.parcel_id.astype(str) != 'nan'
    wd_df.loc[selectedID, 'missinglot'] = wd_df[selectedID].parcel_id.apply(lambda x: without_lots(x))
    return wd_df

def read_taxlot(year, mute=True):
    """
    function to read taxlot data
    """
    txfilepath = os.path.join(txpath, 'Taxlots' + str(year) + '.gdb')
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
    input wd_df is the output from combine_wd_tables (read all set001 data without merging)
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

def combine_taxlot(exportID=True):
    """
    combine taxlots from all years
    """
    frames = []
    for year in range(yearstart, yearend):
        tx_dt = read_taxlot(year)
        tx_dt['year'] = str(year)
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
    counties = get_record_dict(setID, wd_df)[0]
    selected_cnty = df.county.isin(counties[1:])
    record_dict = get_record_dict(setID, wd_df)[1]
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
        ndf = tocheck_df.merge(tdf[['ORTaxlot', 'geometry']], on='ORTaxlot')
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
    if export: 
        # the update will overwrite the first output
        ngdf[~ngdf.geometry.isnull()][selcols].to_file(os.path.join(inpath + f'\\{outfolder}\\', f'matched_records_{setID}.shp'), 
                                                  driver='ESRI Shapefile')  
    return ngdf

def report_unmatched(gdf, setID, nm_to_add, mute = False):
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
    return sorted_df
    
def compare_data_report(gdf, setID, nm_to_add, export = False):
    """
    compare the output data with the existing output from the manual process
    input gdf from match_wd_data_with_taxlot
    return 1) missed_match_ID, missed_gdf: what is missing in the existing matched data; 
           2) addedID, added_gdf: what is added in the existing matched data;
           3) missed_ID: unmatched IDs
    """
    unmatched_wd_df = report_unmatched(gdf, setID, nm_to_add, mute = True)
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
    unmatched_wd_df = report_unmatched(gdf, setID, nm_to_add, mute = True)
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
