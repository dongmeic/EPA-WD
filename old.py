def check_corrected_data(df, setID, all_taxlot, nm_to_add, export=False):
    corrected = pd.read_excel(os.path.join(inpath, 'DSL data originals', 'Data corrections feedback to DSL', 
                           f'DSL Database corrections {setID}.xlsx'))
    corrected.rename(columns={'trsqq': 'cor_trsqq',
                              'parcel_id':'cor_parcel_id'}, inplace=True)
    res = review_with_lots(df, setID, all_taxlot, nm_to_add)
    gdf = res[0]
    df_wlots = res[1]
    IDs1 = check_duplicates(df_wlots.wetdet_delin_number.values)[0]
    IDs2 = check_duplicates(corrected.wetdet_delin_number.values)[0]
    comIDs = [ID for ID in IDs1 if ID in IDs2]
    wdID_to_check = [wdID for wdID in df_wlots.wetdet_delin_number.values if wdID not in corrected.wetdet_delin_number.values]
    cols1 = df_wlots.columns
    cols2 = corrected.columns
    comcols = [col for col in cols1 if col in cols2]
    cols = [col for col in df_wlots.columns if col not in comcols]
    cols.append('wetdet_delin_number')
    df_wlots = df_wlots[cols]
    if len(comIDs) > 0:
        df_wlots = df_wlots[~df_wlots.wetdet_delin_number.isin(comIDs)]
        corrected = corrected[~corrected.wetdet_delin_number.isin(comIDs)]
    cor_df = corrected.merge(df_wlots, on = 'wetdet_delin_number')
    cor_df = cor_df.drop(['trsqq', 'parcel_id'], axis=1)
    cor_df.rename(columns={'cor_trsqq': 'trsqq', 
                          'cor_parcel_id': 'parcel_id',
                          'year': 'Year'}, inplace=True)
    df_wlots_to_check = df_wlots[df_wlots.wetdet_delin_number.isin(wdID_to_check+comIDs)]
    cor_df.loc[:, 'county'] = cor_df.county.apply(lambda x: string.capwords(x))
    cor_df = reindex_data(cor_df)
    cor_df_re = match_wd_data_with_taxlot(df=cor_df, setID=setID, all_taxlot=all_taxlot, nm_to_add=nm_to_add)
    collist = list(gdf.columns)
    collist.remove('recordID')
    ndf = gdf.append(cor_df_re[collist])
    if export:
        ngdf = gpd.GeoDataFrame(ndf, crs="EPSG:2992", geometry='geometry')
        selcols = ['wdID', 'trsqq', 'parcel_id', 'notes', 'lots', 'lot', 'ORTaxlot', 'record_ID', 'geometry']
        ngdf[~ngdf.geometry.isnull()][selcols].to_file(os.path.join(inpath + '\\output\\matched\\', f'matched_records_{setID}.shp'), driver='ESRI Shapefile')
    return ndf, cor_df, comIDs, df_wlots_to_check