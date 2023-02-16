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

inpath = r'L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping'
wdpath = inpath + '\\DSL data originals'
txpath = inpath + '\\GIS\ORMAP_data\ORMAP_Taxlot_Years'
yearstart = 2017
yearend = 2023
#outfolder = 'test'
outfolder = 'output\\matched'
# create a spreadsheet to create a dictionary for the match between county name and code
cnt_ID = pd.read_excel(r'T:\DCProjects\EPA-WD\CNT_Code.xlsx')
# create a dictionary to look up county code
cnt_dict = dict(zip(cnt_ID.COUNTY, cnt_ID.ID))

################################################ Tier 2 #####################################################
def get_point_from_lonlat(lon, lat):
    df = pd.DataFrame([[lon, lat]], columns=['Longitude', 'Latitude'])
    gdf = gpd.GeoDataFrame(df, crs="EPSG:4326", geometry=gpd.points_from_xy(df.Longitude, df.Latitude))
    gdf = gdf.to_crs(epsg=2992)
    return gdf

# point in polygon - WD point in taxtlot
# require taxtlot
def extract_taxlot_info(wd_pt, taxlot, year):
    pip_mask = taxlot.contains(wd_pt.loc[0, 'geometry']) 
    pip_data = taxlot.loc[pip_mask].copy()
    pip_data.loc[:, 'YDiff'] = abs(pip_data.year.astype(int) - int(year))
    ID = pip_data.loc[pip_data.YDiff == np.min(pip_data.YDiff.values), 'ORTaxlot'].values[0]
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

def compare_trsqq(trsqq_to_check, trsqq_to_compare):
    numbers1, letters1 = separate_numbers_letters(trsqq_to_check[:-2])
    numbers2, letters2 = separate_numbers_letters(trsqq_to_compare[:-2])
    letters1.append(trsqq_to_check[-2:])
    letters2.append(trsqq_to_compare[-2:])
    trsqq_to_check_lst = remove_tuple_format(combine_lists(numbers1, letters1))
    trsqq_to_compare_lst = remove_tuple_format(combine_lists(numbers2, letters2))
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

def report_trsqq_correction(trsqq_to_check, trsqq_to_compare):
    diff_idx, correct_trsqq_elements, errors = compare_trsqq(trsqq_to_check, trsqq_to_compare)
    if len(diff_idx) == 1:
        res = f'{trsqq_correction_dict[diff_idx[0]]} is corrected from {errors} to {correct_trsqq_elements}'
    else:
        keylist = [trsqq_correction_dict.get(key) for key in diff_idx]
        if len(keylist) == 2:
            joined_keys = ' and '.join(keylist)
            joined_errors = ' and '.join(errors)
            joined_corrections = ' and '.join(correct_trsqq_elements)
            res = f'{joined_keys} are corrected from {joined_errors} to {joined_corrections} respectively'
        else:
            res = f'{join_list_elements(keylist)} are corrected from {join_list_elements(errors)} to {join_list_elements(correct_trsqq_elements)} respectively'
        return res


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
def get_lot_numbers(text):
    # remove parenthesis from text
    txt = text.replace('(','').replace(')','')
    # split the text
    lot_list = re.split(",|, | |-", txt)
    # remove text elements
    l = []
    # in case there are still number-letter strings (e.g., '1a')
    for t in [lot for lot in lot_list if ~lot.isnumeric()]:
        if any(c.isdigit() for c in t):
            l.append(re.sub('\D', '', t))
    res = unique([lot for lot in lot_list if lot.isnumeric()] + l)
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
    return get_tr_code(x, code='t') + get_tr_code(x, code='r') + get_s_code(x) + '{:0<10}'.format(x)[8] + '{:0<10}'.format(x)[9]

# wd_dt is from read_wd_table or reorganize_tocheck
# return reindexed data
def reindex_data(wd_dt):
    # get a list of lot numbers in each parcel id record
    wd_dt.loc[:, 'lots'] = wd_dt.parcel_id.apply(lambda x: get_lot_numbers(str(x)))
    # repeat the rows based on the number of lot numbers
    ndf = wd_dt.reindex(wd_dt.index.repeat(wd_dt.lots.str.len()))
    # add the column to list the lot for all
    ndf.loc[:, 'lot'] = list(chain.from_iterable(wd_dt.lots.values.tolist()))
    ndf['lots'] = ndf['lots'].apply(lambda x: ', '.join(dict.fromkeys(x).keys()))
    # get county code
    ndf.loc[:, 'cnt_code'] = ndf.county.map(cnt_dict)
    # get OR taxlot IDs for wd data
    ndf.loc[:, 'ORTaxlot'] = ndf[['cnt_code', 'trsqq', 'lot']].apply(lambda row: str(row.cnt_code).zfill(2) + convert_trsqq(row.trsqq) + '--' + ('000000000' + row.lot)[-9:], axis = 1)
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
    wd_dt.loc[selectedID, 'notes'] = wd_dt[selectedID]['parcel_id'].apply(lambda x: make_notes(x))
    ndf = reindex_data(wd_dt)
    # get year from the receive date
    ndf.loc[:, 'recyear'] = ndf.received_date.apply(lambda x: x.year)
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
    wd_df.loc[:, 'record_ID'] = wd_df.loc[:, 'record_ID'] + nm_to_add
    selectedID = wd_df.parcel_id.astype(str) != 'nan'
    wd_df.loc[selectedID, 'notes'] = wd_df[selectedID]['parcel_id'].apply(lambda x: make_notes(x))
    # get year from the receive date
    wd_df.loc[:, 'recyear'] = wd_df.received_date.apply(lambda x: x.year)
    # get year from the wd ID 
    wd_df.loc[:, 'IDyear'] = wd_df.wetdet_delin_number.apply(lambda x: x[2:6]) 
    selectedID = wd_df.parcel_id.astype(str) != 'nan'
    wd_df.loc[selectedID, 'missinglot'] = wd_df[selectedID].parcel_id.apply(lambda x: without_lots(x))
    return wd_df

# function to read taxlot data
def read_taxlot(year, mute=True):
    txfilepath = os.path.join(txpath, 'Taxlots' + str(year) + '.gdb')
    start = time.time()
    tx_dt = gpd.read_file(txfilepath, layer='TL_Dissolv')
    end = time.time()
    if not mute:
        print(f'got taxlot data in {year} and it took about {str(round((end - start)/60, 0))} minutes')
    return tx_dt

# merge wd table with taxlot polygons in a single year
# use this function when checking data by year
def merge_data_by_year(setID, file, year):
    wd_dt = clean_wd_table(setID, file)
    tx_dt = read_taxlot(year)
    merged = wd_dt[wd_dt.IDyear == str(year)].merge(tx_dt[['ORTaxlot', 'geometry']], 
                                                    on='ORTaxlot', 
                                                    how='left')
    #print('got merged data between wd and taxlot')
    return merged

# get a county and record in a dictionary to update recordID
# wd_df is the output from combine_wd_tables (read all set001 data without merging)
def get_record_dict(setID, wd_df):
    counties = wd_df.county.unique()
    count_records = []
    for cnty in counties:
        count = len(wd_df[wd_df.county == cnty])
        count_records.append(count)
    record_df = pd.DataFrame({'county': counties, 'rcount':count_records})
    record_df['cum_count'] = record_df[['rcount']].cumsum(axis = 0, skipna = True).rcount.values
    record_dict = dict(zip(record_df.county[1:len(counties)], record_df.cum_count[0:(len(counties)-1)]))
    return counties, record_dict

# combine taxlots from all years
def combine_taxlot(exportID=False):
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

# update record_ID if needed
# wd_df is the output from combine_wd_tables (read all files in the same set without merging with taxlots)
# df is the reindexed wd data, from combined_reindexed_data
def update_recordID(df, wd_df, setID, nm_to_add):
    counties = get_record_dict(setID, wd_df)[0]
    selected_cnty = df.county.isin(counties[1:])
    record_dict = get_record_dict(setID, wd_df)[1]
    df.loc[selected_cnty, 'record_ID'] = df.loc[selected_cnty, 'recordID'] + df[selected_cnty].county.map(record_dict) + nm_to_add
    df.loc[df.county == counties[0], 'record_ID'] = df.loc[df.county == counties[0], 'recordID'] + nm_to_add
    df.loc[:, 'record_ID'] = df.record_ID.astype('int64', copy=False)
    return df

# combine reindexed wd data in the same set folder
def combined_reindexed_data(setID, nm_to_add):
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

# get the geometry from adjacent years with the same taxlot ID
# df is reindexed from combined_reindexed_data
# return geodata with matched geometry
# run this to include only the matched records with the original ID
def match_wd_data_with_taxlot(df, setID, all_taxlot, nm_to_add, export=False):
    with open(os.path.join(inpath, "ORTaxlot.pkl"), "rb") as f:
        all_txid = pickle.load(f)
    tocheck_txid = df.ORTaxlot.unique()
    found = [txid for txid in tocheck_txid if txid in all_txid]
    if len(found) > 0:
        tocheck_df = df[df.ORTaxlot.isin(found)]
        taxlot_tocheck = all_taxlot[all_taxlot.ORTaxlot.isin(found)]
        taxlot_tocheck = taxlot_tocheck.merge(tocheck_df, on='ORTaxlot', how='left')
        taxlot_tocheck.loc[:, 'ydiff'] = taxlot_tocheck[['year', 'IDyear']].apply(lambda row: abs(int(row.year) - int(row.IDyear)), axis=1)
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
    if export: 
        selcols = ['wdID', 'trsqq', 'parcel_id', 'notes', 'lots', 'lot', 'ORTaxlot', 'record_ID', 'geometry']
        ngdf[~ngdf.geometry.isnull()][selcols].to_file(os.path.join(inpath + f'\\{outfolder}\\', f'matched_records_{setID}.shp'), 
                                                  driver='ESRI Shapefile')  
    return ngdf

# report the unmatched data after matching the taxlot data
# gdf from match_wd_data_with_taxlot
def report_unmatched(gdf, setID, nm_to_add, mute = False):
    wd_df = combine_wd_tables(setID, nm_to_add)
    matched_rID = gdf.record_ID.unique()
    unmatched_wd_df = wd_df[~wd_df.record_ID.isin(matched_rID)]
    nc = unmatched_wd_df[unmatched_wd_df.parcel_id.astype(str) == 'nan'].shape[0]
    nr = round((nc/wd_df.shape[0]) * 100, 2)
    r = round((unmatched_wd_df.shape[0]/wd_df.shape[0]) * 100, 2)
    if not mute:
        print(f'it is about {r}% of data in the original {wd_df.shape[0]} records unmatched')
        print(f'there are {nc} records ({nr}% of the original records) without parcel id')
    return unmatched_wd_df
    
# compare the output data with the existing output from the manual process
# gdf from match_wd_data_with_taxlot
# missed_match_ID, missed_gdf: what is missing in the existing matched data
def compare_data_report(gdf, setID, nm_to_add, export = False):
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

# get the unmatched records with/without taxlot IDs
# gdf from match_wd_data_with_taxlot
def records_with_lots(gdf, setID, nm_to_add, c='Y'):
    unmatched_wd_df = report_unmatched(gdf, setID, nm_to_add, mute = True)
    double_check = unmatched_wd_df[unmatched_wd_df.missinglot == c]
    return double_check

# review the unmatched records with taxlot IDs to rematch with corrected data info
# df is reindexed from combined_reindexed_data
def review_with_lots(df, setID, all_taxlot, nm_to_add):
    n_gdf = match_wd_data_with_taxlot(df, setID, all_taxlot, nm_to_add)
    wd_df = combine_wd_tables(setID, nm_to_add)
    df_wo_lots = records_with_lots(n_gdf, setID, nm_to_add, c='Y')
    wo_lots_ID = df_wo_lots.record_ID.unique()
    to_map_rID = [rID for rID in wd_df.record_ID.values if rID not in n_gdf.record_ID.unique()]
    unmatched_wlots_ID = [rID for rID in to_map_rID if rID not in wo_lots_ID]
    wlots_df = wd_df[wd_df.record_ID.isin(unmatched_wlots_ID)]  
    return n_gdf, wlots_df, to_map_rID

# export mapped data (and missed matched data) from the existing to review the patterns of trsqq
# gdf from match_wd_data_with_taxlot
def review_mapped_data(gdf, setID, all_taxlot, export=False):
    setgdf = gpd.read_file(os.path.join(inpath, 'GIS', 'Join_Statewide.gdb'), layer=f'WD_{setID}_Combined')
    to_map_rID = review_with_lots(setID, all_taxlot)[2]
    mapped = setgdf[setgdf.Record_ID.isin(to_map_rID)]
    if export:
        mapped.rename(columns={'wetdet_delin_number': 'wdID', 
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
                      }, inplace=True)
        selcols = list(mapped.columns[list(map(lambda x: x <= 10, list(map(len, list(mapped.columns)))))])
        mapped[selcols].to_file(os.path.join(inpath + f'\\{outfolder}\\', f'mapped_in_{set_ID}.shp'), 
                                                  driver='ESRI Shapefile')
    return mapped

# update match on the corrected data
# df is reindexed from combined_reindexed_data
# gdf from match_wd_data_with_taxlot
# return matched data with corrected info
def check_corrected_data(df, setID, all_taxlot, nm_to_add, export=False):
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
    cor_df_re = match_wd_data_with_taxlot(df=cor_df, setID=setID, all_taxlot=all_taxlot, nm_to_add=nm_to_add)
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
    with open(os.path.join(inpath, "ORTaxlot.pkl"), "rb") as f:
        taxlotIDs_to_search = pickle.load(f)
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

def read_trsqq():
    with open(os.path.join(inpath, "trsqq_list.pkl"), "rb") as f:
        trsqq = pickle.load(f)
    with open(os.path.join(inpath, "trsqq_dict.pkl"), "rb") as f:
        trsqq_dict = pickle.load(f) 
    df = pd.read_csv(os.path.join(inpath, "trsqq_df.csv"))
    return trsqq, trsqq_dict, df

# get the maybe taxlot
# this is a test function
def get_maybe_taxlot(trsqq_to_check):
    trsqq, trsqq_dict, df = read_trsqq()
    closematch = difflib.get_close_matches(trsqq_to_check, trsqq)
    trsqq_matched = unique(closematch)
    checktaxlot = [*map(trsqq_dict.get, trsqq_matched)]
    string_to_search = trsqq_matched[0][0:8]
    search_res = unique([i for i in unique(df.trsqq) if string_to_search in i])
    if len(checktaxlot) == 1 and len(search_res)==1:
        values = [i for i in trsqq_dict if trsqq_dict[i]==checktaxlot[0]]
    else:
        values = search_res
    return values

# reorganize the tocheck data
# tocheck_df is the last output of check_corrected_data: df_wlots_to_check
# trsqq_n: new trsqq
# n_trsqq: number of possibly corrected trsqq values from the search
# this is a test function
def reorganize_tocheck(tocheck_df):
    tocheck_df.loc[:, 'trsqq_n'] = tocheck_df.loc[:, 'trsqq'].apply(lambda x: get_maybe_taxlot(x))
    tocheck_df.loc[:, 'n_trsqq'] = tocheck_df.loc[:, 'trsqq_n'].apply(lambda x: len(x))
    torematch_df = tocheck_df[tocheck_df.n_trsqq == 1]
    res = get_trsqq_list()
    trsqq_dict = res[1]
    torematch_df.loc[:, 'ORTaxlot'] = torematch_df.trsqq_n.apply(lambda x: [*map(trsqq_dict.get, x)][0])
    return tocheck_df, torematch_df
