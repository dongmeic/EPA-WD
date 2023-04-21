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
1. Run `split_unmatched_df` and `review_unmatched_df_r2`
