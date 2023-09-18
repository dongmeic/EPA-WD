import time
import webbrowser

import pandas as pd

from const import INPATH


class LoopR1Reviewer:
    def __init__(
            self, set_id=None, wd_id_list=None, df=None, is_partial=False,
            do_print_index=False, wd_id=None, wd_df=None):
        self.set_id = set_id
        self.wd_id_list = wd_id_list
        self.df = df
        #self.is_partial = is_partial
        if not is_partial:
            self.df = pd.read_csv(
                fr'{INPATH}\output\to_review\unmatched_df_{set_id}_r1_N.csv')
        if self.df is not None:
            wd_id_list = list(df.wetdet_delin_number.unique())
        else:
            if wd_id_list is None:
                raise ValueError('Noeed to provide a list of wdIDs.')
        self.do_print_index = do_print_index
        self.wd_id = wd_id
        self.wd_df = wd_df

    def review(self):
        'Loop through the unmatched records and check the original records'
        n = len(self.wd_id_list)
        i = -1 if self.wd_id is None else self.wd_id_list.index(self.wd_id)
        to_add, wdID = self._get_wd_ids_to_add(n, i)
        if to_add:
            return to_add, wdID
        else:
            print('Nothing to add in LoopR1Reviewer.review()')

    def _get_wd_ids_to_add(self, n, i):
        to_add = []
        for wdID in self.wd_id_list[i + 1:]:
            print(wdID)
            j = self.wd_id_list.index(wdID)
            remaining = n - j
            print(
                f'{round(((j / n) * 100), 1)}% digitized\n'
                f'{remaining} records remain. '
                f'Estimated time remaining {int(0.8*remaining + 0.5)} hours...')
            if self.do_print_index:
                self._print_index(wdID)
            user_input = input("Press 'p' to pause or any key to stop...")
            if user_input in ['p', 'P']:
                if self._user_validates(wdID):
                    to_add.append(wdID)
            else:
                break
            time.sleep(1)
        return to_add, wdID

    def _print_index(self, wdID):
        if self.df is not None:
            index = (
                self.df[self.df.wetdet_delin_number == wdID].index[0] + 1)
            print(f'index = {index}')
            print(self._check_unmatched_r1(wdID=wdID, df=self.df))
        else:
            print(f'index = {self.wd_id_list.index(wdID) + 1}')
            print(self._check_unmatched_r1(wdID = wdID, df=self.wd_df))

    def _check_unmatched_r1(self, wdID, df):
        'Check unmatched records for a given WDID'
        url = self.df.loc[df.wetdet_delin_number == wdID, 'DecisionLink'].values[0]
        selcols = [
            'county', 'trsqq', 'parcel_id', 'latitude', 'longitude', 'record_ID',
            'notes', 'missinglot', "status_name", "is_batch_file"]
        if str(url) == 'nan':
            print('Decision link is not available')
        else:
            webbrowser.open(url)
        return df.loc[df.wetdet_delin_number == wdID, selcols]

    def _user_validates(wdID):
        while True:
            user_input = input(
                "Press 'a' to add the wd record or 'c' to continue...")
            if user_input in ['a', 'A']:
                return True
            if user_input in ['c', 'C']:
                return False
            
        
