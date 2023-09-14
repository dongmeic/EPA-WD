import os
import re

import fiona
import geopandas as gpd
import pandas as pd

from utils import all_a_in_b, read_geo_data, remove_duplicates


class TaxlotReader:
    def __init__(self, taxlot_path, county_list):
        '''
        Args:
        - taxlot_path (str): path to Taxlot data
        - county_list (list?): list of counties to use
        '''
        self.taxlot_path = taxlot_path
        self.county_list = county_list
        self.list_cols = [
            'County', 'Town', 'TownPart', 'TownDir', 'Range', 'RangePart',
            'RangeDir', 'SecNumber', 'Qtr', 'QtrQtr', 'Anomaly', 'MapSufType',
            'MapNumber', 'ORMapNum', 'Taxlot','MapTaxlot', 'ORTaxlot']
        self.select_cols = [x.capitalize() for x in self.listcols]
        # Next 2 lines not necessary, but help indicate these values will be
        # set as <self> properties in other methods
        self.gdb = None
        self.year = None

    def read(self, year):
        'Read taxlots prior to 2018'
        self.year = year
        self.gdb = self.taxlot_path + f'\\Taxlots{year}.gdb'
        frames = None
        if year in [2016, 2017]:
            frames = self._get_frames_2016_17()
        elif year in [2011, 2014, 2015]:
            frames = self._get_frames_2011_14_15()
        elif year == 2012:
            frames = self._get_frames_2012()
        else:
            raise ValueError(f'Cannot get data for {year}')
        gdf = pd.concat(frames, ignore_index=True)
        gdf['Year'] = year
        return gdf

    def _get_frames_2016_17(self):
        taxlot = gpd.read_file(self.gdb, layer='Taxlots')
        if all([col in taxlot.columns for col in self.list_cols]):
            taxlot = taxlot[self.list_cols]
        else:
            raise ValueError(
                f'check the column names of the {self.year} taxlots!')
        frames = [taxlot]
        return frames

    def _get_frames_2011_14_15(self):
        layer_list = fiona.listlayers(self.gdb)
        selected_layers = [
            x for x in layer_list if x not in list(
                filter(
                    lambda x: re.search(r'_prop_tbl|_TaxCode', x),
                    layer_list))]
        if self.year!= 2011:
            selected_layers = [
                layer for layer in selected_layers
                if layer in self.county_list]
        else:
            selected_layers = [
                layer for layer in selected_layers
                if (layer in self._county_list) and (layer != 'Jefferson')]
        to_remove = [f'{x}_' for x in selected_layers]
        frames = []
        for layer in sorted(selected_layers):
            taxlot = self._get_taxlot_from_layer(layer, to_remove)
            frames.append(taxlot)
        return frames

    def _get_taxlot_from_layer(self, layer, to_remove):
        print('Getting taxlot from layer:', layer)
        taxlot = gpd.read_file(self.gdb, layer=layer)
        taxlot_cols = [re.sub(r'[0-9]', '', x) for x in taxlot.columns]
        col_names = [
            self._remove_county_name(x, to_remove) for x in taxlot_cols]
        taxlot.columns = col_names
        taxlot = (
            self._select_taxlot_col_2011(taxlot, col_names)
            if self.year == 2011
            else self._select_taxlot_cols(taxlot, col_names))
        taxlot = taxlot[self.selected_cols]
        taxlot.columns = self.list_cols
        return taxlot

    def _select_taxlot_cols_2011(self, taxlot, col_names):
        selected_cols = remove_duplicates(col_names)
        if all_a_in_b(['MapTaxlot', 'Maptaxlot'], col_names):
            if 'Maptaxlot' in selected_cols:
                selected_cols.remove('Maptaxlot')
        taxlot = taxlot[selected_cols]
        taxlot.columns = [x.capitalize() for x in selected_cols]
        return taxlot

    def _select_taxlot_cols(self, taxlot, col_names):
        txt = 'taxlot__|taxlot_|Taxlots_|TAXLOT_|Taxlot_|taxlots_|TaxLot_'
        selected_cols = [re.search(txt, x, re.IGNORECASE) for x in col_names]
        selected_cols = list(
            set(selected_cols)
            - {'Taxlots_F', 'Taxlots_Map_Taxlo', 'Taxlots_Accnum'})
        n_col_names = [
            re.sub(txt, '', x, re.IGNORECASE) for x in selected_cols]
        taxlot = taxlot[selected_cols]
        taxlot.columns = [x.capitalize() for x in n_col_names]
        return taxlot

    @staticmethod
    def _remove_county_name(x, to_remove):
        '''Remove county name from strings
        Args:
        - x (list[str]): the input list
        - to_remove (list[str]): list of strings to remove from <x>
        '''
        # Very confusing
        #tf = [tr in x for tr in to_remove]
        #tr = [tr for tr in to_remove if tr in x]
        #if any(tf):
        #    x = re.sub(tr[0], '', x)
        #return x
        x = set(x) - set(to_remove)
        return list(x)

    def _get_frames_2012(self):
        taxlot_dir = self.taxlot_path + f'\\Taxlots{self.year}'
        dir_list = os.listdir(taxlot_dir)
        layers = [
            layer for layer in dir_list if layer not in ['Umatilla', 'Lane']]
        frames = []
        for layer in layers:
            taxlot = self._get_taxlot_from_layer_2012(layer, taxlot_dir)
            frames.append(taxlot)
        return frames

    def _get_taxlot_from_layer_2012(self, layer, taxlot_dir):
        print('Getting taxlots from layer:', layer)
        path = f'{taxlot_dir}\\{layer}'
        txt = 'taxlot|Taxlot'
        files = os.listdir(path)
        file_list = list(
            filter(lambda x: re.search(txt, x, re.IGNORECASE), files))
        selected_files = [f for f in file_list if 'Taxlots' not in f]
        file_path = selected_files[0].split('.')[0]
        taxlot = read_geo_data(f'{path}\\{file_path}.shp')
        taxlot_cols_caps = [x.capitalize() for  x in taxlot.columns]
        if all_a_in_b(self.list_cols, taxlot.columns):
            taxlot = taxlot[self.list_cols]
        elif all_a_in_b(self.select_cols, taxlot_cols_caps):
            taxlot.columns = taxlot_cols_caps
            taxlot = taxlot[self.select_cols]  
            taxlot.columns = self.list_cols
        else:
            print(
                f'Check the column names of the {self.year} taxlots in '
                f'{layer}!')
        return taxlot


class MapIndexReader:
    def __init__(self, taxlot_path, county_dict):
        self.taxlot_path = taxlot_path
        self.county_dict = county_dict
        self.col_names = ['County', 'ORMapNum', 'geometry']
        self.available_years = [2011, 2012]

    def read(self, year):
        if year not in self.available_years:
            raise ValueError(f'Year {year} not available')
        gdf = {2011: self._read_2011, 2012: self._read_2012}[year]
        gdf['Year'] = year
        return gdf

    def _read_2011(self):
        gdb = fr'{self.texlot_path}\Taxlots2011.gdb'
        gdf = gpd.read_file(gdb, layer='mapIndex')
        gdf = gdf[self.col_names]
        return gdf

    def _read_2012(self):
        frames = []
        dir_list = os.listdir(fr'{self.taxlot_path}\Taxlots2012')
        for layer in dir_list:
            layer_gdf = self._get_layer_gdf(layer, 2012)
            frames.append(layer_gdf)
        gdf = pd.concat(frames, ignore_index=True)
        return gdf

    def _get_layer_gdf(self, layer, year):
        path = fr'{self.taxlot_path}\Taxlots{year}\{layer}'
        files = os.listdir(path)
        file_list = list(
            filter(lambda x: re.search('mapindex|MIndex', x, re.IGNORECASE), files))
        file_name = file_list[0].split('.')[0]
        layer_file = fr'{path}\{file_name}.shp'
        gdf = read_geo_data(layer_file)
        if all_a_in_b(self.col_names, gdf.columns):
            gdf = gdf[self.col_names]
        else:
            gdf = self._fix_gdf_columns(gdf, layer, year)
        return gdf

    def _fix_gdf_columns(self, gdf, layer, year):
        gdf.columns = [x.capitalize() for x in gdf.columns]
        column_rename_map = {
            'First_orma': 'ORMapNum', 'Geometry': 'geometry', 'Ormapnum': 'ORMapNum'}
        if all_a_in_b(self.col_names, gdf.columns):
            gdf.rename(columns=column_rename_map, inplace=True)
            gdf = gdf[self.col_names]
        elif 'County' not in gdf.columns:
            gdf['County'] = int(self.county_dict[layer])
            gdf.rename(columns=column_rename_map, inplace=True)
            gdf = gdf[self.col_names]
        else:
            raise ValueError(
                f'Check the column names of the {year} taxmap in '
                f'{layer}!')
        return gdf
