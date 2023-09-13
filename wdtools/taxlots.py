# TODO: continue here--refactor

class TaxlotReader:
    def __init__(self):
        pass

    def read(year):
        'Read taxlots prior to 2017'
        frames = []
        listcols = [
            'County', 'Town', 'TownPart', 'TownDir', 'Range', 'RangePart',
            'RangeDir', 'SecNumber', 'Qtr', 'QtrQtr', 'Anomaly', 'MapSufType',
            'MapNumber', 'ORMapNum', 'Taxlot','MapTaxlot', 'ORTaxlot']
        selcols = list(map(lambda x: x.capitalize(), listcols))
        gdb = txpath + f'\\Taxlots{year}.gdb'
        if year in range(2016, 2018):      
            txlot = gpd.read_file(gdb, layer='Taxlots')
            if all([col in txlot.columns for col in listcols]):
                txlot = txlot[listcols]
            else:
                print(f'check the column names of the {year} taxlots!')
            frames.append(txlot)
        elif year in [2011, 2014, 2015]:
            lyrlist = fiona.listlayers(gdb)
            lyrsel = [
                x for x in lyrlist
                if x not in list(
                    filter(
                        lambda x: re.search(r'_prop_tbl|_TaxCode', x),
                        lyrlist))]
            if year != 2011:
                lyrsel = [lyr for lyr in lyrsel if lyr in cntlst]
            else:
                lyrsel = [
                    lyr for lyr in lyrsel
                    if (lyr in cntlst) and (lyr != 'Jefferson')]
            toRemove = list(map(lambda x: x + '_', lyrsel))
            for lyr in sorted(lyrsel):
                print(lyr)
                txlot = gpd.read_file(gdb, layer=lyr)
                lst = list(
                    map(lambda x: re.sub(r'[0-9]', '', x), txlot.columns))
                colnms = list(map(lambda x: removeCountyNm(x, toRemove), lst))
                txlot.columns = colnms
                if year != 2011:
                    txt = (
                        'taxlot__|taxlot_|Taxlots_|TAXLOT_|Taxlot_|taxlots_|'
                        'TaxLot_')
                    colsel = list(
                        filter(
                            lambda x: re.search(txt, x, re.IGNORECASE), colnms)
                    )
                    colsel = unique(
                        [col for col in colsel if col not in
                         ['Taxlots_F', 'Taxlots_Map_Taxlo', 'Taxlots_Accnum']])
                    ncolnms = list(
                        map(
                            lambda x: re.sub(txt, '', x, re.IGNORECASE),
                            colsel))
                    txlot = txlot[colsel]
                    txlot.columns = list(
                        map(lambda x: x.capitalize(), ncolnms))
                else:
                    if all(
                            [colnm in colnms
                             for colnm in ['MapTaxlot', 'Maptaxlot']]):
                        colsel = [
                            col for col in colnms
                            if col not in (
                                check_duplicates(colnms)[0] + ['Maptaxlot'])]
                    else:
                        colsel = [
                            col for col in colnms
                            if col not in check_duplicates(colnms)[0]]
                    txlot = txlot[colsel]
                    txlot.columns = list(map(lambda x: x.capitalize(), colsel))
                txlot = txlot[selcols]  
                txlot.columns = listcols
                frames.append(txlot)
        elif year == 2012:
            dir_list = os.listdir(txpath + f'\\Taxlots{year}')
            lyrs = [lyr for lyr in dir_list if lyr not in ['Umatilla', 'Lane']]
            for lyr in lyrs:
                print(lyr)
                path = txpath + f'\\Taxlots{year}\\{lyr}'
                txt = 'taxlot|Taxlot'
                files = os.listdir(path)
                filelst = list(
                    filter(lambda x: re.search(txt, x, re.IGNORECASE), files))
                filesel = [file for file in filelst if 'Taxlots' not in file]
                file = filesel[0].split('.')[0]
                txlot = readGeoData(path+f'\\{file}.shp')
                if all([col in txlot.columns for col in listcols]):
                    txlot = txlot[listcols]
                elif all(
                        [col in list(
                            map(lambda x: x.capitalize(), txlot.columns))
                         for col in selcols]):
                    txlot.columns = list(
                        map(lambda x: x.capitalize(), txlot.columns))
                    txlot = txlot[selcols]  
                    txlot.columns = listcols
                else:
                    print(
                        f'check the column names of the {year} taxlots in '
                        f'{lyr}!')
                frames.append(txlot)        
        else:
            print('check the year!')
        gdf = pd.concat(frames, ignore_index=True)
        gdf.loc[:, 'Year'] = year
        return gdf

