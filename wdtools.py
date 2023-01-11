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

# function to clean up wd data
def read_wd_table(setID, file):
    start = time.time()
    datafile = os.path.join(wdpath, setID, file)
    xl = pd.ExcelFile(datafile)
    wd_dt = pd.read_excel(datafile, sheet_name=xl.sheet_names[1])
    wd_dt.loc[:, 'record_ID'] = range(1, wd_dt.shape[0] + 1)
    # get a list of lot numbers in each parcel id record
    wd_dt = wd_dt[wd_dt.parcel_id.astype(str) != 'nan']
    wd_dt.loc[:, 'lots'] = wd_dt.parcel_id.apply(lambda x: get_lot_numbers(x))
    # repeat the rows based on the number of lot numbers
    ndf = wd_dt.reindex(wd_dt.index.repeat(wd_dt.lots.str.len()))
    # add the column to list the lot for all
    ndf.loc[:, 'lot'] = list(chain.from_iterable(wd_dt.lots.values.tolist()))
    # get county code
    ndf.loc[:, 'cnt_code'] = ndf.county.map(cnt_dict)
    # get OR taxlot IDs for wd data
    ndf.loc[:, 'ORTaxlot'] = ndf[['cnt_code', 'trsqq', 'lot']].apply(lambda row: str(row.cnt_code).zfill(2) + row.trsqq[0:2] + '.00' + row.trsqq[2:5] + '.00'+ row.trsqq[5:8] + '{:0<10}'.format(row.trsqq)[8] + '{:0<10}'.format(row.trsqq)[9] + '--' + ('000000000' + row.lot)[-9:], axis = 1)
    # get year from the receive date
    ndf.loc[:, 'year'] = ndf.received_date.apply(lambda x: x.year)
    # get year from the wd ID 
    ndf.loc[:, 'IDyear'] = ndf.wetdet_delin_number.apply(lambda x: x[2:6])
    
    end = time.time()
    #print('cleaned up wd data in {0} and it took about {1} seconds'.format(file, str(end - start)))
    return ndf

def combine_wd_table(setID):
    files = list_files(os.path.join(wdpath, 'Set001'))
    files = [file for file in files if '~$' not in file]
    frames = []
    for file in files:
        datafile = os.path.join(wdpath, setID, file)
        xl = pd.ExcelFile(datafile)
        wd_dt = pd.read_excel(datafile, sheet_name=xl.sheet_names[1])
        frames.append(wd_dt)
    wd_df = pd.concat(frames, ignore_index=True)
    wd_df.loc[:, 'record_ID'] = range(1, wd_df.shape[0] + 1)  
    return wd_df

# function to read taxlot data
def read_taxlot(year):
    txfilepath = os.path.join(txpath, 'Taxlots' + str(year) + '.gdb')
    start = time.time()
    tx_dt = gpd.read_file(txfilepath, layer='TL_Dissolv')
    end = time.time()
    #print('got taxlot data in {0} and it took about {1} minutes'.format(str(year), str(round((end - start)/60, 0))))
    return tx_dt

def merge_data(setID, file, year):
    wd_dt = read_wd_table(setID, file)
    tx_dt = read_taxlot(year)
    merged = wd_dt[wd_dt.IDyear == str(year)].merge(tx_dt[['ORTaxlot', 'geometry']], 
                                                    on='ORTaxlot', 
                                                    how='left')
    #print('got merged data between wd and taxlot')
    return merged

def make_notes(text):
    r = re.search('ROW', text)
    p = re.search('partial', text)
    res = None
    if r and p:
        res = 'ROW & partial'
    elif r:
        res = 'ROW'
    elif p:
        res = 'partial' 
    return res
    