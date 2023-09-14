import json

import fiona
import geopandas as gpd
import pandas as pd
from shapely.geometry import shape


def all_a_in_b(a, b):
    'Determine if all items in iterable a are in iterable b'
    return set(a) - set(b) == set()


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
