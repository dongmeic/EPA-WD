# EPA-WD
Document the tasks and steps on the wetland delineation work funded by Environmental Protection Agency

## Required readings

1. [Oregon Cadastral Data Exchange Standard](https://www.oregon.gov/geo/standards/Cadastral%20Standard%20v3.2.pdf);
2. [Oregon Cadastral Map System](https://digital.osl.state.or.us/islandora/object/osl%3A981082/datastream/OBJ/view)

## Steps

### Tier 1 - initial match with raw data

run `run_Tier1(setID, nm_to_add, all_taxlot)`

Detailed steps:
1. Download and save data in L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping\DSL data originals, by set;
2. Check number to add for the record ID from the previous set, and run `combine_wd_tables`. This will combine all the tables;
3. Check the county for typos and unusual notes, correct typos and leave the unusual records as they are;
4. Run `combined_reindexed_data`, this will reindex all the cleaned records to match with taxlots;
5. Run `match_wd_data_with_taxlot`, this will export the first run of matching;
6. Run `report_unmatched`, this will check how much is matched and export the unmatched records;

### Tier 2 - correct records with parcel ID

1) run `run_Tier2_step1(setID, unmatched_df, all_taxlot)`
2) manual run to get correction notes, `notes_review_setX.ipynb`(X is the set ID)
3) run `run_Tier2_step3(r1_df, r2_df, setID, nm_to_add, wd, all_taxlot)`

Detailed steps:
1. Run `split_unmatched_df` and `review_unmatched_df_r2`;
2. If there is any record reported "the record seems to be correct", review the wdtools functions to get the automatic match;
3. Open "review_unmatched_Set004_r2_N_0.csv" and save it as Excel Workbook, select all and autofit the column width (under format) for better view, this step is only to view the correction comments for notetaking;
4. Mannually create the correction notes on Excel or on Jupter Notebook `notes_review_setX.ipynb` (X is the set ID), which will be saved as "unmatched_df_SetXXX_r2_N_notes.csv" (XXX as the set ID, same as below), can skip the partial taxlots in this step;
5. Review "unmatched_df_SetXXX_r1_N.csv" and create correction notes, can skip the partial taxlots in this step;
6. Run `correct_unmatched`, `combine_corrected_unmatched` (skip `update_unmatched_df_r2`), `reindex_data`, `match_wd_data_with_taxlot`

### Tier 3 & 4 - digitize the partial and unmatched records

Digitize the unmatched and partially matched records - `digitize_setX_loop.ipynb`

### Report
1. Combine digitized and matched records, including WD ID;
2. Review unmatched or issue IDs;
3. Get the deliverables including the matched or digitized WD records;
4. Report the issues for final digitizing work in needed.

## Notes

The scripts are organized to run by set. A set is collected when DSL sends data at certain time following their master schedule. 

## WD SA data description

### Objectives and processes

DSL requested data attritutes and test data for their database development tests on June 2023. The process involves discussion on the sample data, steps to review study area polygons to get sample data, and follow-up on the criteria for final deliverables.

### Study area polygons

Certain issues need to be considered: