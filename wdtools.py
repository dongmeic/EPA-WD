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

inpath = r'L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping'
wdpath = inpath + '\\DSL data originals'
txpath = inpath + '\\GIS\ORMAP_data\ORMAP_Taxlot_Years'
yearstart = 2017
yearend = 2023
outfolder = 'test'

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
    wd_dt.loc[:, 'record_ID'] = range(1, wd_dt.shape[0] + 1)
    return wd_dt

def reindex_data(wd_dt):
    # get a list of lot numbers in each parcel id record
    wd_dt.loc[:, 'lots'] = wd_dt.parcel_id.apply(lambda x: get_lot_numbers(x))
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
    wd_dt = wd_dt[wd_dt.parcel_id.astype(str) != 'nan']
    wd_dt.loc[:, 'notes'] = wd_dt.parcel_id.apply(lambda x: make_notes(x))
    ndf = reindex_data(wd_dt)
    # get year from the receive date
    ndf.loc[:, 'year'] = ndf.received_date.apply(lambda x: x.year)
    # get year from the wd ID 
    ndf.loc[:, 'IDyear'] = ndf.wetdet_delin_number.apply(lambda x: x[2:6]) 
    end = time.time()
    #print(f'cleaned up wd data in {file} and it took about {end - start} seconds')
    return ndf

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
def combine_wd_table(setID):
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
    wd_df.loc[:, 'year'] = wd_df.received_date.apply(lambda x: x.year)
    # get year from the wd ID 
    wd_df.loc[:, 'IDyear'] = wd_df.wetdet_delin_number.apply(lambda x: x[2:6]) 
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
def merge_data(setID, file, year):
    wd_dt = clean_wd_table(setID, file)
    tx_dt = read_taxlot(year)
    merged = wd_dt[wd_dt.IDyear == str(year)].merge(tx_dt[['ORTaxlot', 'geometry']], 
                                                    on='ORTaxlot', 
                                                    how='left')  # to include all the records in the original table, using "left"
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

# get a county and record in a dictionary
def get_record_dict(setID):
    # read all set001 data without merging
    wd_df = combine_wd_table(setID)
    selectedID = wd_df.parcel_id.astype(str) != 'nan'
    wd_df.loc[selectedID, 'notes'] = wd_df[selectedID]['parcel_id'].apply(lambda x: make_notes(x))
    counties = wd_df.county.unique()
    count_records = []
    for cnty in counties:
        count = len(wd_df[wd_df.county == cnty])
        count_records.append(count)
    record_df = pd.DataFrame({'county': counties, 'rcount':count_records})
    record_df['cum_count'] = record_df[['rcount']].cumsum(axis = 0, skipna = True).rcount.values
    record_dict = dict(zip(record_df.county[1:len(counties)], record_df.cum_count[0:(len(counties)-1)]))
    return counties, record_dict

# merge data by year and combine all tables in the same set
def combine_merged_data(setID, export=False):
    files = list_files(os.path.join(wdpath, setID))
    files = [file for file in files if '~$' not in file]
    start1 = time.time()
    frames = []
    # need to change the year start and end when they are changed
    for year in range(yearstart, yearend):
        tx_dt = read_taxlot(year)
        for file in files:
            wd_dt = clean_wd_table(setID, file = file)  
            df = wd_dt[wd_dt.IDyear == str(year)].merge(tx_dt[['ORTaxlot', 'geometry']], 
                                                        on='ORTaxlot')
            frames.append(df)
    df = pd.concat(frames, ignore_index=True)      
    end1 = time.time()
    print(f'it took {round((end1 - start1)/60, 0)} minutes to complete {setID}')
    df['response_date'] = df['response_date'].dt.strftime("%Y-%m-%d")
    df['received_date'] = df['received_date'].dt.strftime("%Y-%m-%d")
    df['lots'] = df['lots'].apply(lambda x: ' '.join(dict.fromkeys(x).keys()))
    gdf = gpd.GeoDataFrame(df, crs="EPSG:2992", geometry='geometry')
    counties = get_record_dict(setID)[0]
    selected_cnty = gdf.county.isin(counties[1:])
    record_dict = get_record_dict(setID)[1]
    gdf.loc[selected_cnty, 'nrecordID'] = gdf.loc[selected_cnty, 'record_ID'] + gdf[selected_cnty].county.map(record_dict)
    gdf.loc[gdf.county == counties[0], 'nrecordID'] = gdf.loc[gdf.county == counties[0], 'record_ID']
    gdf.loc[:, 'nrecordID'] = gdf.nrecordID.astype('int64', copy=False)
    if export:
        gdf.rename(columns={'wetdet_delin_number': 'wdID', 
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
        gdf[~gdf.geometry.isnull()].to_file(os.path.join(inpath + f'\\{outfolder}\\', f'matched_records_{setID}.shp'), 
                                                  driver='ESRI Shapefile')  
    return gdf

# report the unmatched data after matching the same year data
def report_unmatched(gdf, setID, mute = False):
    wd_df = combine_wd_table(setID)
    matched_rID = gdf.nrecordID.unique()
    unmatched_wd_df = wd_df[~wd_df.record_ID.isin(matched_rID)]
    nc = unmatched_wd_df[unmatched_wd_df.parcel_id.astype(str) == 'nan'].shape[0]
    nr = round((nc/wd_df.shape[0]) * 100, 2)
    r = round((unmatched_wd_df.shape[0]/wd_df.shape[0]) * 100, 2)
    if not mute:
        print(f'it is about {r}% of data in the original {wd_df.shape[0]} records unmatched')
        print(f'there are {nc} records ({nr}% of the original records) without parcel id')
    return unmatched_wd_df
    
# compare the output data with the existing output from the manual process
def compare_data_report(gdf, setID, export = False):
    unmatched_wd_df = report_unmatched(gdf, setID, mute = True)
    # unmatched IDs from the run
    missed_ID = unmatched_wd_df.record_ID.unique()
    setgdf = gpd.read_file(os.path.join(inpath, 'GIS', 'Join_Statewide.gdb'), layer=f'WD_{setID}_Combined')
    matched_rID = gdf.nrecordID.unique()
    # missed IDs in the existing data that is not nan
    missed_gdf = setgdf[setgdf.Record_ID.astype(str) != 'nan'][~setgdf.Record_ID.isin(matched_rID)]
    missedID = missed_gdf.Record_ID.astype('int64', copy=False)
    # matched IDs in the existing data that is not nan
    matched_gdf = setgdf[(setgdf.Record_ID.astype(str) != 'nan') & (setgdf.Record_ID.isin(matched_rID))]
    matchedID = matched_gdf.Record_ID.unique()
    missed_match_ID = [ID for ID in list(matched_rID) if ID not in list(matchedID)]
    missed_gdf = gdf[gdf.nrecordID.isin(missed_match_ID)]
    if export:
        missed_gdf[~missed_gdf.geometry.isnull()].to_file(os.path.join(inpath + f'\\{outfolder}\\', f'missed_records_in_{setID}_res.shp'), 
                                                  driver='ESRI Shapefile')
    return missed_match_ID, missed_gdf

# combine taxlots
def combine_taxlot():
    frames = []
    for year in range(yearstart, yearend):
        tx_dt = read_taxlot(year)
        tx_dt['year'] = year
        frames.append(tx_dt[['year', 'ORTaxlot', 'geometry']])
    df = pd.concat(frames, ignore_index=True)
    gdf = gpd.GeoDataFrame(df, crs="EPSG:2992", geometry='geometry')
    return gdf

# get the geometry from adjacent years with the same taxlot ID
def rematch_data(gdf, setID, all_taxlot):
    all_txid = all_taxlot.ORTaxlot.unique()
    double_check = records_with_lots(gdf, setID, c='N')
    double_check_ri = reindex_data(double_check)
    tocheck_txid = double_check_ri.ORTaxlot.unique()
    found = [txid for txid in tocheck_txid if txid in all_txid]
    if len(found) > 0:
        tocheck_df = double_check_ri[double_check_ri.ORTaxlot.isin(found)]
        taxlot_tocheck = all_taxlot[all_taxlot.ORTaxlot.isin(found)]
        taxlot_tocheck = taxlot_tocheck.merge(tocheck_df[['record_ID', 'ORTaxlot', 'IDyear']], on='ORTaxlot', how='left')
        taxlot_tocheck.loc[:, 'ydiff'] = taxlot_tocheck[['year', 'IDyear']].apply(lambda row: abs(row.year - int(row.IDyear)), axis=1)
        df = taxlot_tocheck.sort_values(by=['ORTaxlot', 'ydiff'])
        # keep the taxlot from the closest year
        df = df.drop_duplicates(subset='ORTaxlot', keep="first")
        tocheck_gdf = tocheck_df.merge(df[['ORTaxlot', 'geometry']], on='ORTaxlot')
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

# merge the matched from exact taxlot IDs
def merge_matched(setID, all_taxlot, export=False):
    gdf = gpd.read_file(os.path.join(inpath + f'\\{outfolder}\\', f'matched_records_{setID}.shp'), driver='ESRI Shapefile')
    tocheck_gdf = rematch_data(gdf=gdf, setID=setID, all_taxlot=all_taxlot)
    gdf = gdf.drop(['record_ID'], axis=1)
    gdf.rename(columns={'nrecordID': 'record_ID'}, inplace=True)
    tocheck_gdf = tocheck_gdf[gdf.columns]
    n_gdf = gdf.append(tocheck_gdf)
    if export:
        n_gdf['lots'] = n_gdf['lots'].apply(lambda x: ' '.join(dict.fromkeys(x).keys()))
        n_gdf[['record_ID', 'wdID', 'parcel_id', 'ORTaxlot', 'geometry']].to_file(os.path.join(inpath + f'\\{outfolder}\\', f'combined_records_in_{setID}.shp'), driver='ESRI Shapefile')
    return n_gdf

# get the unmatched records with/without taxlot IDs
def records_with_lots(gdf, setID, c='Y'):
    unmatched_wd_df = report_unmatched(gdf, setID, mute = True)
    double_check = unmatched_wd_df[unmatched_wd_df.missinglot == c]
    return double_check

# export mapped data (and missed matched data) from the existing to review the patterns of trsqq
def review_mapped_data(gdf, setID, all_taxlot, export=False):
    setgdf = gpd.read_file(os.path.join(inpath, 'GIS', 'Join_Statewide.gdb'), layer=f'WD_{setID}_Combined')
    to_map_rID = review_with_lots(gdf, setID, all_taxlot)[1]
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

# review the unmatched records with taxlot IDs to rematch with a revised trsqq
def review_with_lots(gdf, setID, all_taxlot):
    n_gdf = merge_matched(setID, all_taxlot)
    wd_df = combine_wd_table(setID)
    df_wo_lots = records_with_lots(gdf, setID, c='Y')
    wo_lots_ID = df_wo_lots.record_ID.unique()
    to_map_rID = [rID for rID in wd_df.record_ID.values if rID not in n_gdf.record_ID.unique()]
    unmatched_wlots_ID = [rID for rID in to_map_rID if rID not in wo_lots_ID]
    wlots_df = wd_df[wd_df.record_ID.isin(unmatched_wlots_ID)]  
    return wlots_df, to_map_rID
