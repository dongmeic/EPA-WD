import re

from utils import remove_duplicates


def get_lot_numbers(parcel_id):
    return LotNumberGetter().get(parcel_id)
    

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
