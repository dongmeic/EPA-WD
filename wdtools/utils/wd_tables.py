class WDTableCombiner:
    def __init__(setID, n_to_addd, is_raw=True):
        '''Combine all the wd tables in the set to review unique records, record_ID
        is used for combined tables.
        Use this function when reindex is not neccessary
        Args:
        -n_to_add (int) the number of previous records (records from the previous
        sets; same for all functions with this variable)
    '''
