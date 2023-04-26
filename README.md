# EPA-WD
Document the tasks and steps on the wetland delineation work funded by Environmental Protection Agency

## Required readings

1. [Oregon Cadastral Data Exchange Standard](https://www.oregon.gov/geo/standards/Cadastral%20Standard%20v3.2.pdf);
2. [Oregon Cadastral Map System](https://digital.osl.state.or.us/islandora/object/osl%3A981082/datastream/OBJ/view)

## Steps

### Tier 1 - initial match with raw data
1. Download and save data in L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping\DSL data originals, by set;
2. Check number to add for the record ID from the previous set, and run `combine_wd_tables`. This will combine all the tables;
3. Check the county for typos and unusual notes, correct typos and leave the unusual records as they are;
4. Run `combined_reindexed_data`, this will reindex all the cleaned records to match with taxlots;
5. Run `match_wd_data_with_taxlot`, this will export the first run of matching;
6. Run `report_unmatched`, this will check how much is matched and export the unmatched records;

### Tier 2 - correct records with parcel ID
1. Run `split_unmatched_df` and `review_unmatched_df_r2`;
2. If there is any record reported "the record seems to be correct", review the wdtools functions to get the automatic match;
3. Open "review_unmatched_Set004_r2_N_0.csv" and save it as Excel Workbook, select all and autofit the column width (under format) for better view;
4. Mannually create the correction notes on Excel or on Jupter Notebook `notes_review_set4`, which will be saved as "unmatched_df_SetXXX_r2_N_notes.csv" (XXX as the set ID), can skip the partial taxlots in this step;
5. 
