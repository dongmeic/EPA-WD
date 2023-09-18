import os
import pickle
import re

import geopandas as gpd
import pandas as pd

from const import (
    COL_DICT, COUNTY_DICT, INPATH, OR_COUNTIES, TID_DST_0, TRANSFORMER, WD_PATH)
from utils.reindexing import reindex_data
from utils.utils import convert_trsqq, list_files


def combine_wd_tables(setID, n_to_add, is_raw=True):
    return WDTableHandler(setID).combine(n_to_add, is_raw)


def clean_wd_table(setID, file_name):
    return WDTableHandler(setID).clean(file_name)
                   

class WDTableHandler:
    def __init__(self, setID):
        self.setID = setID

    def combine(self, n_to_add, is_raw):
        '''Combine all the wd tables in the set to review unique records, record_ID
        is used for combined tables.
        Use this function when reindex is not neccessary
        Args:
        -n_to_add (int) the number of previous records (records from the previous
        sets; same for all functions with this variable)
        '''
        if is_raw:
            wd_df = self._get_wd_df_from_raw()
        else:
            wd_df = pd.read_csv(fr'{WD_PATH}\Corrected_by_Set\{self.setID}.csv')
        # this creates unique IDs for all the records in the same set
        wd_df['record_ID'] = range(1, wd_df.shape[0] + 1)
        wd_df['record_ID'] = wd_df.copy().record_ID + n_to_add
        selected_ids = wd_df.parcel_id.astype(str) != 'nan'
        wd_df.loc[selected_ids, 'notes'] = (
            wd_df[selected_ids]['parcel_id'].apply(lambda x: self._make_notes(x)))
        # get year from the receive date
        if self.is_raw:
            wd_df.loc[:, 'recyear'] = wd_df.received_date.apply(lambda x: x.year)
        else:
            wd_df.loc[:, 'recyear'] = wd_df.received_date.apply(
                lambda x: int(x.split('-')[0]))
        # get year from the wd ID
        wd_df['IDyear'] = wd_df.wetdet_delin_number.apply(lambda x: x[2:6])
        selected_ids = wd_df.parcel_id.astype(str) != 'nan'
        wd_df.loc[selected_ids, 'missinglot'] = (
            wd_df[selected_ids].parcel_id.apply(
                lambda x: self._is_without_lots(x)))
        return wd_df

    def _get_wd_df_from_raw(self):
        frames = []
        files = list_files(os.path.join(WD_PATH, self.setID))
        # in case there are unidentified files
        files = [f for f in files if '~$' not in f]
        for f in files:
            data_file = os.path.join(WD_PATH, self.setID, f)
            xl = pd.ExcelFile(data_file)
            wd_dt = pd.read_excel(data_file, sheet_name=xl.sheet_names[1])
            frames.append(wd_dt)
        wd_df = pd.concat(frames, ignore_index=True)
        return wd_df

    def clean(self, file_name):
        'Clean up wd data from a single file'
        wd_dt = self._read_wd_table(file_name)
        # this will help identify problematic records with numbers
        selected_IDs = wd_dt.parcel_id.astype(str) != 'nan'
        wd_dt.loc[selected_IDs, 'notes'] = (
            wd_dt.copy()[selected_IDs]['parcel_id'].apply(
                lambda x: self._make_notes(x)))
        wd_dt['county'] = wd_dt['county'].apply(lambda x: x.title())
        wd_dt = wd_dt[wd_dt.county.isin(OR_COUNTIES)]
        if wd_dt.empty:
            return wd_dt
        else:
            wd_dt = self._clean_nonempty_wd_table(wd_dt, selected_IDs)
            return wd_dt

    def _clean_nonempty_wd_table(self, wd_dt, selected_IDs):
        ndf = reindex_data(wd_dt)
	# get year from the receive date
        ndf.loc[:, 'recyear'] = ndf.copy().received_date.apply(
            lambda x: x.year)
        # get year from the wd ID
        ndf.loc[:, 'IDyear'] = ndf.wetdet_delin_number.apply(lambda x: x[2:6])
        ndf.loc[selected_IDs, 'missinglot'] = (
            ndf.loc[selected_IDs, 'parcel_id'].apply(
                lambda x: self._is_without_lots(x)))
        ndf['response_date'] = ndf['response_date'].dt.date
        ndf['received_date'] = ndf['received_date'].dt.date
        ndf['county'] = ndf['county'].apply(lambda x: x.title())
        res = ndf[ndf.county.isin(OR_COUNTIES)]
        return res

    def _read_wd_table(self, file_name):
        'Read wd tables; recordID is used for single tables'
        data_file = os.path.join(WD_PATH, self.setID, file_name)
        xl = pd.ExcelFile(data_file)
        wd_dt = pd.read_excel(data_file, sheet_name=xl.sheet_names[1])
        wd_dt.loc[:, 'county'] = wd_dt.county.apply(lambda x: x.capitalize())
        # this will show the records in the original single table and the ID will
        # be up dated when the tables are combined       
        wd_dt['recordID'] = range(1, wd_dt.shape[0] + 1)
        return wd_dt

    @staticmethod
    def _make_notes(text):
        'Add notes to the wd records'
        r = re.search('ROW|RR', text, re.IGNORECASE)
        p = re.search('partial|part|p|portion', text, re.IGNORECASE)
        m = re.search('Many|multiple|SEVERAL|various', text, re.IGNORECASE)
        w = re.search('Water', text, re.IGNORECASE)
        l = re.search('RAIL', text, re.IGNORECASE)
        res = []
        for srch, msg in zip(
                [r, p, m, w, l], ['ROW', 'Partial', 'Many', 'Water', 'Rail']):
            if srch:
                res.append(msg)
        res = ', '.join(res)
        return res

    @staticmethod
    def _is_without_lots(text):
        '''Check whether the parcel ID is without lots
        Returns (str): 'Y'/'N'
        '''
        if any(c.isdigit() for c in text):
            names = re.findall(r'\d+', text)
            suffixes = ['st', 'nd', 'rd', 'th']
            matches = (
                [f'{name}{suff}' for name in names for suff in suffixes]
                + [f'{name} {suff}' for name in names for suff in suffixes])
            has_matches = all([any([x in text for x in matches])])
            return 'Y' if has_matches else 'N'
        else:
            return 'Y'


class WDSAJoiner:
    def __init__(
            self, df, gdf, map_index, out_name='wd_mapped_data', do_export=True):
        '''
        Args:
        - df (pandas.DataFrame): corrected WD dataframe from get_all_wd()
        - gdf (geopanda.GeoDataFrame): the SA polygons from get_all_SA()
        - map_index (type?): the combined taxmaps from read_all_mapIdx()
        '''
        self.df = df
        self.gdf = gdf
        self.map_index = map_index
        self.out_name = out_name
        self.do_export = do_export

    def join(self):
        self.df['ORMapNum'] = self.df[['county', 'trsqq']].apply(
            lambda row: self._create_ORMap_name(
                county=row.county, trsqq=row.trsqq),
            axis = 1)
        frames = self._get_frames_from_wd_list()
        sa_df = pd.concat(frames, ignore_index=True)
        sa_gdf = gpd.GeoDataFrame(sa_df, geometry='geometry')
        if self.do_export:
            self._export(sa_gdf)
        return sa_gdf
    
    @staticmethod
    def _create_ORMap_name(county, trsqq):
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

    def _get_frames_from_wd_list(self):
        wd_list = self.gdf.wdID.unique()
        frames = []
        for wid in wd_list:
            df_wid = self.df[self.df.wetdet_delin_number == wid]
            gdf_wid = self.gdf[self.gdf.wdID == wid]
            or_map_n = df_wid.ORMapNum.values
            yr = wid[2:6]
            map_idx = self.map_index[
                (self.map_index.ORMapNum.isin(or_map_n))
                & (self.map_index.year == yr)]
            exp_gdf = self._split_WD_to_taxmaps(
                df=df_wid, gdf=gdf_wid, wdID=wid, mapindex=map_idx)
            exp_gdf['lat'], exp_gdf['lon'] = TRANSFORMER.transform(
                exp_gdf.representative_point().x, exp_gdf.representative_point().y)
            ndf = df_wid.drop_duplicates(subset='ORMapNum')
            ndf.drop(columns=['parcel_id','site_id','record_ID'], inplace=True)
            g = df_wid.groupby('ORMapNum')
            pi_df, si_df, ri_df = [
                g
                .apply(lambda x: '; '.join(x[col].astype(str).unique()))
                .to_frame(name='parcel_id')
                .reset_index(drop=True)
                for col in ['parcel_id', 'site_id', 'record_id']]
            sdf = pd.concat([ri_df, pi_df, si_df], axis=1)
            fdf = ndf.merge(sdf, on='ORMapNum')
            exp_gdf = fdf.merge(
                exp_gdf[['code', 'ORMapNum', 'lat', 'lon', 'geometry']],
                on='ORMapNum')
            frames.append(exp_gdf)
        return frames

    @staticmethod
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

    def _export(self, sa_gdf):
        sa_gdf = sa_gdf.rename(columns=COL_DICT)
        path =  os.path.join(INPATH, 'output', 'final', f'{self.out_name}.shp')
        try:
            sa_gdf.to_file(path, index=False)
        except RuntimeError:
            sa_gdf['geometry'] = sa_gdf.geometry.buffer(0)
            sa_gdf.to_file(path, index=False)
        print(f'Wrote sa_gdf to {path}')
