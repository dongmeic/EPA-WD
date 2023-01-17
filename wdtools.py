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

inpath = r'L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping'
wdpath = inpath + '\\DSL data originals'
txpath = inpath + '\\GIS\ORMAP_data\ORMAP_Taxlot_Years'
yearstart = 2017
yearend = 2023
outfolder = 'test'

# gdf below generally refers to the matched records
# create a spreadsheet to create a dictionary
cnt_ID = pd.read_excel(r'T:\DCProjects\EPA-WD\CNT_Code.xlsx')

# create a dictionary to look up county code
cnt_dict = dict(zip(cnt_ID.COUNTY, cnt_ID.ID))

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

# read wd tables, record_ID is used for single tables
def read_wd_table(setID, file):
    datafile = os.path.join(wdpath, setID, file)
    xl = pd.ExcelFile(datafile)
    wd_dt = pd.read_excel(datafile, sheet_name=xl.sheet_names[1])
    # this will show the records in the original single table and the ID will be updated when the tables are combined
    wd_dt.loc[:, 'recordID'] = range(1, wd_dt.shape[0] + 1)
    return wd_dt

def reindex_data(wd_dt):
    # get a list of lot numbers in each parcel id record
    wd_dt.loc[:, 'lots'] = wd_dt.parcel_id.apply(lambda x: get_lot_numbers(str(x)))
    # repeat the rows based on the number of lot numbers
    ndf = wd_dt.reindex(wd_dt.index.repeat(wd_dt.lots.str.len()))
    # add the column to list the lot for all
    ndf.loc[:, 'lot'] = list(chain.from_iterable(wd_dt.lots.values.tolist()))
    # get county code
    ndf.loc[:, 'cnt_code'] = ndf.county.map(cnt_dict)
    # get OR taxlot IDs for wd data
    ndf.loc[:, 'ORTaxlot'] = ndf[['cnt_code', 'trsqq', 'lot']].apply(lambda row: str(row.cnt_code).zfill(2) + row.trsqq[0:2] + '.00' + row.trsqq[2:5] + '.00'+ row.trsqq[5:8] + '{:0<10}'.format(row.trsqq)[8] + '{:0<10}'.format(row.trsqq)[9] + '--' + ('000000000' + row.lot)[-9:], axis = 1)
    return ndf

# function to clean up wd data
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
    ndf['lots'] = ndf['lots'].apply(lambda x: ' '.join(dict.fromkeys(x).keys()))
    end = time.time()
    #print(f'cleaned up wd data in {file} and it took about {end - start} seconds')
    return ndf

# check whether the parcel ID that is without lots
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
def combine_wd_tables(setID):
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

# add notes to the wd records
def make_notes(text):
    r = re.search('ROW', text)
    p = re.search('partial|part|p', text)
    m = re.search('Many|MANY|multiple|many', text)
    res = None
    if r and p and m:
        res = 'ROW & partial & many'
    elif r and p:
        res = 'ROW & partial'
    elif r and m:
        res = 'ROW & many'
    elif p and m:
        res = 'partial & many'
    elif r:
        res = 'ROW'
    elif p:
        res = 'partial'
    elif m:
        res = 'many'    
    return res

# get a county and record in a dictionary to update recordID
def get_record_dict(setID):
    # read all set001 data without merging
    wd_df = combine_wd_tables(setID)
    selectedID = wd_df.parcel_id.astype(str) != 'nan'
    wd_df.loc[selectedID, 'notes'] = wd_df.loc[selectedID, 'parcel_id'].apply(lambda x: make_notes(x))
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
        with open(os.path.join(inpath, 'test', "ORTaxlot.pkl"), "wb") as f:
            pickle.dump(list(gdf.ORTaxlot.unique()), f) 
    return gdf

# update record_ID if needed
def update_recordID(df, setID):
    counties = get_record_dict(setID)[0]
    selected_cnty = df.county.isin(counties[1:])
    record_dict = get_record_dict(setID)[1]
    df.loc[selected_cnty, 'record_ID'] = df.loc[selected_cnty, 'recordID'] + df[selected_cnty].county.map(record_dict)
    df.loc[df.county == counties[0], 'record_ID'] = df.loc[df.county == counties[0], 'recordID']
    df.loc[:, 'record_ID'] = df.record_ID.astype('int64', copy=False)
    return df

# get the geometry from adjacent years with the same taxlot ID
# df is reindexed from combined_reindexed_data
# return geodata with matched geometry
def match_wd_data_with_taxlot(df, setID, all_taxlot):
    all_txid = all_taxlot.ORTaxlot.unique()
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
        tocheck_gdf = tocheck_df.merge(tdf[['ORTaxlot', 'geometry']], on='ORTaxlot')
        tocheck_gdf.rename(columns={'wetdet_delin_number': 'wdID', 
                      'address_location_desc':'loc_desc', 
                      'Coord-Source': 'coord_src',
                      'DocumentName':'doc_name',
                      'DecisionLink':'doc_link',
                      'is_batch_file':'isbatfile',
                      'status_name': 'status_nm',
                      'received_date':'receiveddt', 
                      'response_date':'responsedt',
                      'reissuance_response_date':'reissuance', 
                      }, inplace=True)
        return tocheck_gdf

# combine reindexed wd data
def combined_reindexed_data(setID):
    frames = []
    files = list_files(os.path.join(wdpath, setID))
    # in case there are unidentified files
    files = [file for file in files if '~$' not in file]
    for file in files:
        wd_dt = clean_wd_table(setID, file = file)  
        frames.append(wd_dt)
    df = pd.concat(frames, ignore_index=True)
    wd_df = update_recordID(df, setID)
    return wd_df 
    
# combine all tables in the same set with the match in the same or adjacent year
def match_taxlot(setID, all_taxlot, export=False):
    all_txid = all_taxlot.ORTaxlot.unique()
    wd_df = combined_reindexed_data(setID)
    gdf = match_wd_data_with_taxlot(df=wd_df, setID=setID, all_taxlot=all_taxlot)
    if export:
        gdf[~gdf.geometry.isnull()].to_file(os.path.join(inpath + f'\\{outfolder}\\', f'matched_records_{setID}.shp'), 
                                                  driver='ESRI Shapefile')  
    return gdf

# report the unmatched data after matching the same year data
# gdf from match_taxlot
def report_unmatched(gdf, setID, mute = False):
    wd_df = combine_wd_tables(setID)
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
# gdf from match_taxlot
def compare_data_report(gdf, setID, export = False):
    unmatched_wd_df = report_unmatched(gdf, setID, mute = True)
    # unmatched IDs from the run
    missed_ID = unmatched_wd_df.record_ID.unique()
    setgdf = gpd.read_file(os.path.join(inpath, 'GIS', 'Join_Statewide.gdb'), layer=f'WD_{setID}_Combined')
    matched_rID = gdf.record_ID.unique()
    # missed IDs in the existing data that is not nan
    missed_gdf = setgdf[setgdf.Record_ID.astype(str) != 'nan'][~setgdf.Record_ID.isin(matched_rID)]
    missedID = missed_gdf.Record_ID.astype('int64', copy=False)
    # matched IDs in the existing data that is not nan
    matched_gdf = setgdf[(setgdf.Record_ID.astype(str) != 'nan') & (setgdf.Record_ID.isin(matched_rID))]
    matchedID = matched_gdf.Record_ID.unique()
    missed_match_ID = [ID for ID in list(matched_rID) if ID not in list(matchedID)]
    missed_gdf = gdf[gdf.record_ID.isin(missed_match_ID)]
    if export:
        missed_gdf[~missed_gdf.geometry.isnull()].to_file(os.path.join(inpath + f'\\{outfolder}\\', f'missed_records_in_{setID}_res.shp'), 
                                                  driver='ESRI Shapefile')
    return missed_match_ID, missed_gdf

# get the unmatched records with/without taxlot IDs
# gdf from match_taxlot
def records_with_lots(gdf, setID, c='Y'):
    unmatched_wd_df = report_unmatched(gdf, setID, mute = True)
    double_check = unmatched_wd_df[unmatched_wd_df.missinglot == c]
    return double_check

# review the unmatched records with taxlot IDs to rematch with corrected data info
# gdf from match_taxlot
def review_with_lots(setID, all_taxlot):
    n_gdf = match_taxlot(setID, all_taxlot)
    wd_df = combine_wd_tables(setID)
    df_wo_lots = records_with_lots(n_gdf, setID, c='Y')
    wo_lots_ID = df_wo_lots.record_ID.unique()
    to_map_rID = [rID for rID in wd_df.record_ID.values if rID not in n_gdf.record_ID.unique()]
    unmatched_wlots_ID = [rID for rID in to_map_rID if rID not in wo_lots_ID]
    wlots_df = wd_df[wd_df.record_ID.isin(unmatched_wlots_ID)]  
    return n_gdf, wlots_df, to_map_rID

# export mapped data (and missed matched data) from the existing to review the patterns of trsqq
# gdf from match_taxlot
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
# gdf from match_taxlot
# return matched data with corrected info
def check_corrected_data(setID, all_taxlot, export=False):
    corrected = pd.read_excel(os.path.join(inpath, 'DSL data originals', 'Data corrections feedback to DSL', 
                           f'DSL Database corrections {setID}.xlsx'))
    corrected.rename(columns={'trsqq': 'cor_trsqq',
                              'parcel_id':'cor_parcel_id'}, inplace=True)
    res = review_with_lots(setID, all_taxlot)
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
    ngdf = gdf.append(cor_df_re[collist])
    if export:
        ngdf.to_file(os.path.join(inpath + f'\\{outfolder}\\', f'combined_records_in_{setID}.shp'), driver='ESRI Shapefile')
    return ngdf, cor_df, comIDs, df_wlots_to_check

# update match with similar trsqq
def get_trsqq_list():
    with open(os.path.join(inpath, 'test', "ORTaxlot.pkl"), "rb") as f:
        taxlotIDs_to_search = pickle.load(f)
    taxlotIDs_cleaned = unique([txID for txID in taxlotIDs_to_search if len(txID) == 29])
    trsqq = list(map(lambda x: x[2:4] + x[7:10] + x[13:18], taxlotIDs_cleaned))
    trsqq_dict = dict(zip(trsqq, taxlotIDs_cleaned))
    df = pd.DataFrame({'ORTaxlot':taxlotIDs_cleaned, 'trsqq':trsqq})
    return trsqq, trsqq_dict, df

# get the maybe taxlot
def get_maybe_taxlot(trsqq_to_check):
    res = get_trsqq_list()
    trsqq = res[0] 
    trsqq_dict = res[1]
    df = res[2]
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
# tocheck_df is the last output of check_corrected_data - df_wlots_to_check
# trsqq_n - new trsqq
# n_trsqq - number of possibly corrected trsqq values from the search
def reorganize_tocheck(tocheck_df):
    tocheck_df.loc[:, 'trsqq_n'] = tocheck_df.loc[:, 'trsqq'].apply(lambda x: get_maybe_taxlot(x))
    tocheck_df.loc[:, 'n_trsqq'] = tocheck_df.loc[:, 'trsqq_n'].apply(lambda x: len(x))
    return tocheck_df
