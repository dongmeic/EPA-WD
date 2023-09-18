from itertools import chain
import re

from ..const import ALL_TXID, COUNTY_DICT, TID_DST, TID_DST_0, TID_DST_1
from lot_numbers import get_lot_numbers
from trsqq_conversion import convert_trsqq


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
