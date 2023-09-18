import os
import re

import pandas as pd

from ..const import OR_COUNTIES, WD_PATH
from reindexing import reindex_data
from util import list_files


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
