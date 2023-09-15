from collections import Counter
from itertools import chain
import json
import os
import pickle
from random import sample
import re

import fiona
import geopandas as gpd
import pandas as pd
from shapely.geometry import shape

from const import ALL_TXID, COUNTY_DICT, INPATH, TID_DST, TID_DST_0, TID_DST_1


def all_a_in_b(a, b):
    'Determine if all items in iterable a are in iterable b'
    return set(a) - set(b) == set()


def check_duplicates(v):
    '''Check for duplicates in a list
    Args:
    - v: (list)
    Returns:
    - (list, int): list of duplicate items, number of items duplicated         
    '''
    item_list = [
        item for item, count in Counter(v).items() if count > 1]
    return item_list, len(item_list)


def convert_trsqq(x):
    '''Convert township, range, section, and quarter-quarter to the taxlot id
    format'''
    TrsqqConverter().convert(x)


class TrsqqConverter:
    def __init__(self):
        self.directions = list('NSEW')

    def convert(self, x):
        x = f'{x:<08s}'
        xt = (
            self._get_tr_code(x)
            + self._get_tr_code(x)
            + self._get_s_code(x)
            + self._get_qq_code(x))
        return xt[:16]

    def _get_tr_code(self, x):
        'Get township/range code from the taxlot'
        nms = re.findall('\d+', x)
        nm1 = nms[1][0:2]
        lts = re.findall("[a-zA-Z]+", x)
        if lts[1] in self.directions:
            lts2 = lts[1]
        else:
            lts2 = sample(self.directions, 1)[0]
        if len(nm1) == 1:
            tr1 = f'0{nm1}'
        elif len(nm1) == 3:
            tr1 = nm1[1:3]
        else:
            tr1 = nm1
        if ('V' in lts2) or ('Y' in lts2):
            tr2 = '.50'
            tr3 = lts2[1]
        elif ('X' in lts2) or ('Z' in lts2):
            if any([x in lts2 for x in ['XS', 'ZN', 'XE', 'ZW']]):
                tr2 = '.75'
            else:
                tr2 = '.25'
        else:
            tr2 = '.00'
            tr3 = lts2
        res = tr1 + tr2 + tr3
        return res

    @staticmethod
    def _get_s_code(x):
        'Get section code from trsqq code'
        if len(x) <= 7:
            s = '00'
        else:
            nms = re.findall('\d+', x)
            n = len(nms)
            k = len(nms[1])
            if (n < 3) and (k > 2):
                nm1 = nms[1][(k-2):(k+1)]
            else:
                nm1 = nms[2]
            if len(nm1) == 1:
                s = f'0{nm1}'
            else:
                s = nm1
        return s

    def _get_qq_code(self, x):
        'Get QQ code from trsqq code'
        nms = re.findall('\d+', x)
        lts = re.findall("[a-zA-Z]+", x)
        if (len(lts) == 2) and (lts[1] not in self.directions):
            if len(lts[1]) == 2:
                qq = lts[1]
            else:
                qq = f'{lts[1]}0'
        elif len(nms[2]) > 2:
            if nms[2] == 4:
                qq = nms[2][2:4]
            else:
                qq = f'{nms[2][2]}0'
        elif len(lts) == 3:
            if len(lts[2]) == 2:
                qq = lts[2]
            else:
                qq = f'{lts[2]}0'
        else:
            if len(nms[0]) == 1:
                t = f'0{nms[0]}'
            else:
                t = nms[0]
            if len(nms[1]) == 1:
                r = f'0{nms[1]}'
            else:
                r = nms[1]
            if len(nms[2]) == 1:
                s = f'0{nms[2]}'
            else:
                s = nms[2]
            trsqq = t + lts[0] + r + lts[1] + s
            qq = f'{trsqq:0<10}'[8:10]
        return qq


def create_ORMap_name(county, trsqq):
    'Return ORMap number based on county name and trsqq'
    with open(os.path.join(INPATH, 'ORMapIndex.pkl'), 'rb') as f:
        all_map_idx = pickle.load(f)
    part1 = str(int(COUNTY_DICT[county])).zfill(2) + convert_trsqq(trsqq)
    map_idx = f'{part1}--0000'
    if map_idx not in all_map_idx:
        if part1 in TID_DST_0:
            for mid in [part1+f'--{x}000' for x in ['D', 'S', 'T']]:
                if mid in all_map_idx:
                    map_idx = mid
    return map_idx


def get_lot_numbers(parcel_id):
    LotNumberGetter().get(parcel_id)


class LotNumberGetter:
    def get(self, parcel_id):
        'Takes a string or integer and returns a list of lot numbers'
        if parcel_id is None:
            return None
        if type(parcel_id) is int:
            res = self._get_from_int(parcel_id)
        else:
            res = self._get_from_str(parcel_id)
        return res

    @staticmethod
    def _get_from_int(parcel_id):
        s = str(parcel_id)
        if len(str(parcel_id)) > 5:
            idx = [i for i, char in enumerate(s) if char != '0']
            lot_list = []
            for i in range(len(idx) - 1):
                lot_list.append(s[idx[i]:idx[i + 1]])
            lot_list.append(s[idx[len(idx) - 1]:])
            res = lot_list
        else:
            res = [s]
        return res

    def _get_from_str(self, parcel_id):
        # remove parenthesis from text 
        if '(' in str(parcel_id):
            txt = parcel_id.replace('(', '').replace(')', '')
        else:
            txt = str(parcel_id)
        lot_list = self._get_lot_list_from_str(txt)
        numbers = self._remove_text_elements(lot_list)
        res = remove_duplicates(
            [lot for lot in lot_list if lot.isnumeric()] + numbers)
        road = re.search('ROW|RR', parcel_id, re.IGNORECASE)
        water = re.search('Water', parcel_id, re.IGNORECASE)
        rail = re.search('RAIL', parcel_id, re.IGNORECASE)
        for srch, msg in zip([road, water, rail], ['ROADS','WATER','RAILS']):
            if srch:
                res.append(msg)
        return res

    @staticmethod
    def _get_lot_list_from_str(txt):
        # split the text
        lot_list = []
        for r in re.split(",|, | ", txt):
            if '-' in r:
                start, end = r.split('-')
                if start.isdigit() and end.isdigit():
                    lot_list += list(range(int(start), int(end) + 1))
                    lot_list = [str(x) for x in lot_list]
            else:
                lot_list.append(r)
        return lot_list

    @staticmethod
    def _remove_text_elements(lot_list):
        numbers = []
        # in case there are still number-letter strings (e.g., '1a')
        for t in [lot for lot in lot_list if ~lot.isnumeric()]:
            if any(c.isdigit() for c in t):
                numbers.append(re.sub('\D', '', t))
        return numbers


def read_geo_data(layer_file):
    '''Read geodata with a geometry check
    Args:
    - layer_file (str): path to layer file
    '''
    try:
        gdf = gpd.read_file(layer_file)
    except ValueError:
        collection = list(fiona.open(layer_file, 'r'))
        df = pd.DataFrame(collection)
        
        # Check Geometry
        def is_valid(geom):
            try:
                shape(geom)
                return 1
            except:
                return 0
            
        df['isvalid'] = df['geometry'].apply(lambda x: is_valid(x))
        df = df[df['isvalid'] == 1]
        collection = json.loads(df.to_json(orient='records'))
        # Convert to geodataframe
        gdf = gpd.GeoDataFrame.from_features(collection)
    return gdf


def reindex_data(wd_dt):
    return DataReindexer.reindex(wd_dt)


class DataReindexer:
    def reindex(self, wd_dt):
        '''Reindex the data based on the number of lots in each parcel id    
        Args:
        - wd_dt (type?):  from read_wd_table or reorganize_tocheck
        make sure the records are with lot numbers  
        Returns (type?) reindexed data
        '''
        # get a list of lot numbers in each parcel id record
        # will need to review the records without any parcel ids
        selected_id = wd_dt.parcel_id.astype(str) != 'nan'
        wd_dt = wd_dt.copy()[selected_id]
        wd_dt.loc[:, 'lots'] = wd_dt['parcel_id'].apply(
            lambda x: get_lot_numbers(x))
        # repeat the rows based on the number of lot numbers
        ndf = wd_dt.reindex(wd_dt.index.repeat(wd_dt.lots.str.len()))
        # add the column to list the lot for all
        ndf.loc[:, 'lot'] = list(
            chain.from_iterable(wd_dt.lots.values.tolist()))
        ndf.loc[:, 'lots'] = ndf['lots'].apply(
            lambda x: ', '.join(dict.fromkeys(x).keys()))
        # get county code
        ndf.loc[:, 'cnt_code'] = ndf.county.map(COUNTY_DICT)
        # get OR taxlot IDs for wd data
        ndf = ndf[~ndf.cnt_code.isnull()]
        ndf.loc[:, 'ORTaxlot'] = ndf[['cnt_code', 'trsqq', 'lot']].apply(
            lambda row: self._create_OR_taxlot(
                cnt_code=row.cnt_code, trsqq=row.trsqq, lot=row.lot),
            axis = 1)
        return ndf

    def _create_OR_taxlot(self, cnt_code, trsqq, lot):
        '''Create the taxlot id based on the county code, township, range,
        section, and lot number
        '''
        part1 = str(int(cnt_code)).zfill(2) + convert_trsqq(trsqq)
        taxlot_id = f'{part1}--{("000000000" + lot)[-9:]}'
        tid_dst_2 = [x[-len(lot):] for x in TID_DST_1]
        if (taxlot_id not in ALL_TXID):
            taxlot_id = self._get_taxlot_id(part1, lot, tid_dst_2)
        return taxlot_id

    def _get_taxlot_id(part1, lot, tid_dst_2):
        if (part1 in TID_DST_0) and (lot in tid_dst_2):
            taxlot_id_list = [
                tid for tid in TID_DST
                if (re.search(part1, tid)) and (re.search(lot, tid))]
            if len(taxlot_id_list):
                taxlot_id = taxlot_id_list[0]
            else:
                n_taxlot_id = (
                    f'{taxlot_id[:6]}{taxlot_id[7:12]}{taxlot_id[13:18]}'
                    f'00--{taxlot_id.split("--")[1]}')
                if n_taxlot_id in ALL_TXID:
                    taxlot_id = n_taxlot_id
        return taxlot_id


def remove_duplicates(lst):
    'Remove duplicates from a list <lst>'
    return list(set(lst))


def split_WD_to_taxmaps(df, gdf, wdID, map_index):
    '''Splits the WD SA ploygons by taxmap                           
    Args:
    - gdf (GeoDatFrame): geodataframe that contains the selected WD ID
    - map_index (GeoDataFrame): taxmap geodataframe of the year
    '''
    selected_map_ids = df[df.wetdet_delin_number==wdID].ORMapNum.unique()
    gdf = gdf[gdf.wdID == wdID]
    selected_map_index = map_index[['ORMapNum','geometry']][
        map_index.ORMapNum.isin(selected_map_ids)]
    try:
        inter = gpd.overlay(
            gdf, selected_map_index, how='intersection', keep_geom_type=False)
    except NotImplementedError:
        gdf['geometry'] = gdf.geometry.buffer(0)
        inter = gpd.overlay(
            gdf, selected_map_index, how='intersection', keep_geom_type=False)
    inter = inter.dissolve('ORMapNum')
    inter['ORMapNum'] = inter.index
    inter.reset_index(drop=True, inplace=True)
    return inter
