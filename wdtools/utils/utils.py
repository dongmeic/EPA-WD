from collections import Counter
import json
import os
import pickle

import fiona
import geopandas as gpd
import pandas as pd
from shapely.geometry import shape

from const import COUNTY_DICT, INPATH, TID_DST_0
from trsqq_conversion import TrsqqConverter


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


def list_files(path, is_folder=False):
    '''This function takes a path and returns a list of files in the path and its
    subdirectories
    '''
    f = []
    for (dirpath, dirnames, filenames) in os.walk(path):
        f.extend(filenames)
        break
    if is_folder:
        f = [x[0] for x in os.walk(path)]
    return f

    
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
