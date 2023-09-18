import os
import pickle
import re

import pandas as pd
from pyproj import Transformer


INPATH = r'L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping'
OUTPATH = fr'{INPATH}\output'
TAXLOT_PATH = fr'{INPATH}\GIS\ORMAP_data\ORMAP_Taxlot_Years'
WD_PATH = fr'{INPATH}\DSL data originals'

# create a spreadsheet to create a dictionary for the match between county
# name and code
COUNTY_IDS = pd.read_excel(fr'{INPATH}\notes\CNT_Code.xlsx')
# create a dictionary to look up county code
COUNTY_DICT = dict(zip(COUNTY_IDS.COUNTY, COUNTY_IDS.ID))
OR_COUNTIES = list(COUNTY_DICT.keys())
SELECTED_VARS = [
    'wetdet_delin_number', 'trsqq', 'parcel_id','address_location_desc',
    'city', 'county', 'site_name', 'site_desc','latitude', 'longitude',
    'DocumentName', 'DecisionLink','is_batch_file', 'status_name',
    'received_date', 'response_date','reissuance_response_date', 'project_id',
    'site_id', 'SetID','record_ID', 'ORMapNum']
VAR_LIST = SELECTED_VARS + ['geometry', 'code']
COL_DICT = {
    'address_location_desc': 'loc_desc',
    'Coord-Source': 'CordSource',
    'DecisionLink': 'doc_link',
    'DocumentName': 'doc_name',
    'is_batch_file': 'isbatfile',
    'received_date': 'receiveddt',
    'reissuance_response_date': 'reissuance',
    'response_date': 'responsedt',
    'status_name': 'status_nm',
    'wetdet_delin_number': 'wdID'}
TRANSFORMER = Transformer.from_crs('EPSG:2992', 'EPSG:4326')
YEAR_START, YEAR_END = 2016, 2023


def read_trsqq():
    'Read trsqq list, dictionary, and dataframe'
    with open(os.path.join(INPATH, 'trsqq_list.pkl'), 'rb') as f:
        trsqq = pickle.load(f)
    with open(os.path.join(INPATH, 'trsqq_dict.pkl'), 'rb') as f:
        trsqq_dict = pickle.load(f)
    df = pd.read_csv(os.path.join(INPATH, 'trsqq_df.csv'))
    return trsqq, trsqq_dict, df


TRSQQ, TRSQQ_DICT, TTDF = read_trsqq()
TID_DST = [
    tid for tid in TTDF.ORTaxlot.unique()
    if any(substring in tid for substring in ['--D', '--S', '--T'])]
TID_DST_0  = [re.split('--', x)[0] for x in TID_DST]
TID_DST_1 = [re.split('--', x)[1] for x in TID_DST]
TSQ_DST = TTDF[TTDF.ORTaxlot.isin(TID_DST)].trsqq.unique()

with open(os.path.join(INPATH, 'ORTaxlot.pkl'), 'rb') as f:
    ALL_TXID = pickle.load(f)

    
