import geopandas as gpd
import pandas as pd

from const import (
    COL_DICT, OR_COUNTIES, OUTPATH, SELECTED_VARS, TRANSFORMER, VAR_LIST)
from utils import (
    all_a_in_b, check_duplicates, create_ORMap_name, reindex_data,
    split_WD_to_taxmaps)


class SASplitter:
    def __init__(
            self, wd_df, sa_gdf_all, all_map_idx, all_taxlot, em_wids,
            do_export=False, out_name='example_data', do_review=False):
        ''' Split study area polygons where multiple record IDs exist.
        rid is record ID.
        Args:
        - wd_df (pandas.DataFrame?): dataframe that includes the example WD IDs
        - sa_gdf_all (GeoDataFrame?): geodataframe that includes the example WD
          IDs
        - all_map_idx (type?): the combined geodataframe of map index from all
          yearsn
        - all_taxlot (GeoDataFrame?): the combined geodataframe of taxlots from
          all years
        - em_wids (type?): the example WD IDs
        - do_export (bool): export if true
        - out_name (str): name of output file
        - do_review (bool): check if all county names are valid
        Returns:
        - (GeoDataFrame?): the combined geodataframe from split_WD_to_records
        '''
        self.wd_df = wd_df
        self.sa_gdf_all = sa_gdf_all
        self.all_map_idx = all_map_idx
        self.all_taxlot = all_taxlot
        self.em_wids = em_wids
        self.do_export = do_export
        self.out_name = out_name
        self.do_review = do_review

    def split_by_rid(self):
        wd_df_s = self.wd_df[self.wd_df.wetdet_delin_number.isin(self.em_wids)]
        self.wd_df_s['ORMapNum'] = wd_df_s[['county', 'trsqq']].apply(
            lambda row: create_ORMap_name(county=row.county, trsqq=row.trsqq),
            axis = 1)
        self.sa_gdf_s = self.sa_gdf_all[
            self.sa_gdf_all.wdID.isin(self.em_wids)]
        self.wdID_list, n = check_duplicates(
            self.wd_df_s.wetdet_delin_number.values)
        if n:
            gdf = self._get_gdf_from_frames()
        else:
            gdf = self._get_gdf_from_merge()
        gdf['lat'], gdf['lon'] = TRANSFORMER.transform(
            gdf.representative_point().x, gdf.representative_point().y)
        if self.do_export:
            self._export(gdf)
        return gdf

    def _get_gdf_from_frames(self):
        frames = []
        for wid in self.wdID_list:
            if wid in self.sa_gdf_s.wdID.values:
                frame = self._get_frame_from_wid(wid)
            frames.append(frame)
        df = pd.concat(frames, ignore_index=True)
        gdf1 = gpd.GeoDataFrame(df, crs="EPSG:2992", geometry='geometry')
        wdID_sdf = self.wd_df_s[
            ~self.wd_df_s.wetdet_delin_number.isin(self.wdID_list)]
        wdID_delin_list = wdID_sdf.wetdet_delin_number.values
        sa_sgdf = self.sa_gdf_s[self.sa_gdf_s.wdID.isin(wdID_delin_list)]
        sa_sgdf.rename(columns={'wdID': 'wetdet_delin_number'}, inplace=True)
        gdf2 = wdID_sdf.merge(
            sa_sgdf[['code', 'wetdet_delin_number', 'geometry']],
            on='wetdet_delin_number')
        gdf = pd.concat([gdf1[VAR_LIST], gdf2[VAR_LIST]])
        return gdf
    
    def _get_frame_from_wid(self, wid):
        or_map_n = (
            self.wd_df_s[self.wd_df_s.wetdet_delin_number == wid]
            .ORMapNum
            .values)
        yr = wid[2:6]
        map_idx = self.all_map_idx[
            (self.all_map_idx.ORMapNum.isin(or_map_n))
            & (self.all_map_idx.year == yr)]
        taxlot = self.all_taxlot[self.all_taxlot.year == yr]
        frame = self._split_WD_to_records(
            df=self.wd_df_s,
            gdf=self.sa_gdf_s,
            wdID=wid,
            mapindex=map_idx,
            taxlots=taxlot,
            do_review=self.do_review)
        return frame

    def _split_WD_to_records(
            self, df, gdf, wdID, map_index, taxlots, do_review=False):
        '''Splits the WD SA ploygons to ploygons by records
        Args:
        - df (pandas.DataFrame): dataframe containing the selected WD ID and
          wetdet_delin_number in the columns
        - gdf (GeoDataFrame?): geodataframe containing the selected WD ID
        - map_index (GeoDataFrame?): taxmap geodataframe of the year
        - taxlots (type?) are the taxlots of the year
        - do_review (bool): check if all county names are valid
        '''
        ndf = df[df.wetdet_delin_number == wdID]
        counties = ndf.county.unique()
        if not self._counties_are_valid(counties, wdID):
            return
        if 'ORMapNum' not in ndf.columns:
            ndf['ORMapNum'] = ndf[['county', 'trsqq']].apply(
                lambda row: create_ORMap_name(
                    county=row.county, trsqq=row.trsqq),
                axis=1)
        trsqq_list, n = check_duplicates(ndf.trsqq.values)
        ngdf = split_WD_to_taxmaps(
            df=ndf, gdf=gdf, wdID=wdID, mapindex=map_index)
        if n:
            out = self._get_df_from_trsqq_list(trsqq_list, ndf, ngdf, taxlots)
        else:
            out = ndf.merge(
                ngdf[['ORMapNum', 'geometry', 'code']],
                on='ORMapNum',
                how='left')
            out = out[VAR_LIST]
        out = gpd.GeoDataFrame(out, geometry='geometry')
        out = out.dissolve('record_ID')
        out['record_ID'] = out.index
        out.reset_index(drop=True, inplace=True)
        return out

    def _counties_are_valid(self, counties, wdID):
        if (len(counties) > 1) and not self.do_review:
            print(f'{wdID} crosses counties!')
            return False
        ## ?? I think this is a syntax err
        #elif counties.any() not in OR_COUNTIES:
        elif not all_a_in_b(counties, OR_COUNTIES):
            print(f'Check the county name {counties[0]}!')
            return False
        return True

    def _get_df_from_trsqq_list(self, trsqq_list, ndf, ngdf, taxlots):
        frames = []
        for trsqq in trsqq_list:
            inter = self._split_taxmap_to_records(
                df=ndf, gdf=ngdf, trsqq=trsqq, taxlots=taxlots)
            frames.append(inter)
        # idf is intersection data frame
        idf = pd.concat(frames, ignore_index=True)
        idf = idf[['geometry', 'record_ID']].merge(
            ndf[SELECTED_VARS],on='record_ID', how='left')
        idf['code'] = ngdf.code.values[0]
        # odf is the other dataframe that excludes intersection taxmaps
        ogdf = ngdf[~ngdf.ORMapNum.isin(idf.ORMapNum)][
            ['ORMapNum', 'geometry', 'code']]
        odf = ndf[~ndf.ORMapNum.isin(idf.ORMapNum)][SELECTED_VARS]
        odf = ogdf.merge(odf, on='ORMapNum', how='left')
        out = pd.concat([idf[VAR_LIST], odf[VAR_LIST]])
        return out

    def _split_taxmap_to_records(df, gdf, trsqq, taxlots):
        '''Applied when the same taxmap has multiple records
        Args:
        - df (pandas.DataFrame): dataframe of the selected WD ID
        - gdf (GeoDataFrame): geodataframe of the selected taxmaps, from
          split_WD_to_taxmaps
        - trsqq (type?): the selected trsqq that appears in multiple record IDs
        - taxlots (type?):  the taxlots of the year
        '''
        df = df[df.trsqq == trsqq]
        rdf = reindex_data(df)
        taxlots = taxlots[taxlots.ORTaxlot.isin(rdf.ORTaxlot.unique())]
        t_df = rdf[['ORTaxlot','record_ID']].merge(
            taxlots[['ORTaxlot', 'geometry']], on='ORTaxlot', how='left')
        t_gdf = gpd.GeoDataFrame(t_df, geometry='geometry')
        t_gdf = t_gdf.dissolve('record_ID')
        t_gdf['record_ID'] = t_gdf.index
        i_gdf = gdf[gdf.ORMapNum==df.ORMapNum.unique()[0]]
        i_gdf = i_gdf.dissolve('wdID')
        i_gdf['wdID'] = i_gdf.index
        i_gdf.reset_index(drop=True, inplace=True)
        inter = gpd.overlay(
            i_gdf, t_gdf, how='intersection', keep_geom_type=False)
        return inter

    def _get_gdf_from_merge(self):
        self.sa_gdf_s.rename(
            columns={'wdID': 'wetdet_delin_number'}, inplace=True)
        gdf = self.wd_df_s.merge(
            self.sa_gdf_s[['code', 'wetdet_delin_number', 'geometry']],
            on='wetdet_delin_number')
        gdf = gdf[VAR_LIST]
        gdf = gpd.GeoDataFrame(gdf, crs="EPSG:2992", geometry='geometry')
        return gdf

    def _export(self, gdf):
        gdf = gdf.rename(columns=COL_DICT)
        path = fr'{OUTPATH}\test\{self.out_name}.shp'
        try:
            gdf.to_file(path)
        except RuntimeError:
            gdf['geometry'] = gdf.geometry.buffer(0)
            gdf.to_file(path)
        print('Saved gdf to', path)
