# EPA-WD
Document the tasks and steps on the wetland delineation work funded by Environmental Protection Agency

## Required readings

1. [Oregon Cadastral Data Exchange Standard](https://www.oregon.gov/geo/standards/Cadastral%20Standard%20v3.2.pdf);
2. [Oregon Cadastral Map System](https://digital.osl.state.or.us/islandora/object/osl%3A981082/datastream/OBJ/view)

## ORMAP Taxlots

The ORMAP taxlots from 2016 to 2023 are ready to use. However, the ORMAP taxlots before 2016 require data cleaning to keep the schema consistent. Data from some years in some counties miss the key information such as ORTaxlot, MapTaxlot, or ORMapNum and Taxlot, and is considerated unusuable. Saving data in geodatabase is easier for better ArcGIS Pro performance. 

## Steps

### Steps on the scripts

Whenever there is "set#" in the script name, make a copy of the script with the current set number to keep the records. This is also because that each set may have different data structure that requires updates on the functions, which is a debugging process. The Jupyter Notebooks are not cleaned scripts and require careful review if they will be used.

1. Run the script `run_and_review_set#.ipynb`

Start with combing tables and scanning the TRSQQ (correct the errors first if necessary), then run the function `run_Tier1` to get the initial matches. Run the function `run_Tier2_step1` to get the notes for the corrections in the next step. Pause and move to the next step.

2. Run the script `notes_review_set#.ipynb`

This step basically creates the notes to correct the WD records based on coordinates. Need to make sure that the coordinates are correct and the study area matched from the taxlots is consistent with the study area from the decision link. 

3. Continue the script `run_and_review_set#.ipynb`

Run the function `run_Tier2_step3` to get the additional matches. Skip this step if the step `run_Tier2_step1` is skipped.

4. Correct the WD data based on the errata

Once the errata is collected (saved at 'output\corrected\corrected_Set#.xlsx' based on 'unmatched_df_Set#_r1_N_notes.csv' and 'unmatched_df_Set#_r2_N_notes.csv' from `notes_review_set#.ipynb`), correct the WD tables after this or next step using the script [`08_correct_WD_records`](https://github.com/dongmeic/EPA-WD/blob/main/08_correct_WD_records.ipynb) based on the errata and save an updated copy ('DSL data originals\Corrected_by_Set\Set#.csv'). It is considered as Tier2 Step2. The corrected WD data will be used in the QAQC step.

5. Run the script `QAQC_run_the_loop.ipynb` for QAQC

The QAQC step reviews the matched SA polygons and the WD decision link on the study area, and check whether they are matched and organize questions for DSL if the decision link is unclear. 

6. Run the script `digitize_set#_loop.ipynb`

First, review the matched records without notes on the parcel ids such as partial taxlots and roads using the function `review_loop_r1` to identify the QAQC records. Then update the [application](https://lcog.maps.arcgis.com/apps/instant/charts/index.html?appid=69fe51df1ce544e980e27e5a5a89dd06) with the QAQC records, matched records without notes, unmatched records, and issue IDs. The data can be updated by replacing the [feature layer](https://lcog.maps.arcgis.com/home/item.html?id=2a9bcd28a8e34516b9f91f312864d544). The next step is digitizing the unmatched or partially-matched records. It is considered as Tier 3 and 4 with feedback on errata and issue IDs for DSL. 

When start digitizing, create a file geodatabase 'L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping\GIS\ArcGIS Pro Project\DataReview\Set*.gdb' and two feature classes 'Set*_wo_lot' and 'Set*_partial'. The notes 'Set*_edited.txt' and 'Set*_edited_1.txt' in the folder 'L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping\output\matched' are created to document the edited taxlots in the shapefile 'matched_records_Set*_edited.shp' or edited taxlots in the file geodatabase. Make a copy of the output file 'matched_records_Set#.shp' created from `run_and_review_set#.ipynb` and rename it as 'matched_records_Set#_edited.shp' to edit the WD SA polygons that are partially matched with taxlots and only require minor editing. 

7. Run the script [`04_combine_matched_digitized`](https://github.com/dongmeic/EPA-WD/blob/main/04_combine_matched_digitized.ipynb)

This step applies `run_Tier3_4_final` and export the first version of final WD SA data.

8. Run the script [`09_add_issue_IDs_to_set`](https://github.com/dongmeic/EPA-WD/blob/012078acc20779062fbe918a5974ba3c11775ceb/09_add_issue_IDs_to_set.ipynb)

This step finalizes the WD SA data with issue IDs. 

9. Update the application data using `update_QAQC_data_on_the_application.ipynb`[https://github.com/dongmeic/EPA-WD/blob/main/update_QAQC_data_on_the_application.ipynb]

The original table was prepared using the script [organize_data_for_QAQC_reporting.ipynb](https://github.com/dongmeic/EPA-WD/blob/main/organize_data_for_QAQC_reporting.ipynb) and updated using the script [update_QAQC_data_on_the_application.ipynb](https://github.com/dongmeic/EPA-WD/blob/main/update_QAQC_data_on_the_application.ipynb). Basically, we need to update the feature layer "[Wetland Delineation and Determination Counts](https://lcog.maps.arcgis.com/home/item.html?id=2a9bcd28a8e34516b9f91f312864d544)" (source data - L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping\reporting\WD_Counts.zip). The function `update_QAQC_data` is created to update the source data for the feature layer. The [web map](https://lcog.maps.arcgis.com/apps/mapviewer/index.html?webmap=c3bffa7cd230464db051ec263b5e4517) and [application](https://lcog.maps.arcgis.com/apps/instant/charts/index.html?appid=69fe51df1ce544e980e27e5a5a89dd06) are both automatically updated when the feature layer is updated.

### Steps on the functions
#### Tier 1 - initial match with raw data

run `run_Tier1(setID, nm_to_add, all_taxlot)`

Detailed steps:
1. Download and save data in L:\NaturalResources\Wetlands\Local Wetland Inventory\WAPO\EPA_2022_Tasks\Task 1 WD Mapping\DSL data originals, by set;
2. Check number to add for the record ID from the previous set, and run `combine_wd_tables`. This will combine all the tables;
3. Check the county for typos and unusual notes, correct typos and leave the unusual records as they are;
4. Run `combined_reindexed_data`, this will reindex all the cleaned records to match with taxlots;
5. Run `match_wd_data_with_taxlot`, this will export the first run of matching;
6. Run `report_unmatched`, this will check how much is matched and export the unmatched records;

#### Tier 2 - correct records with parcel ID

This step can be skipped if there is not unmatched data with specific parcel IDs. 

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

#### Tier 3 & 4 - digitize the partial and unmatched records

Digitize the unmatched and partially matched records - `digitize_setX_loop.ipynb`

1. Run `run_Tier3_4_final`(`04_combine_matched_digitized.ipynb`) to combine digitized and matched records, including WD ID;
2. Review unmatched or issue IDs;
3. Get the deliverables including the matched or digitized WD records;
4. Report the issues to DSL for the final digitizing work needed;
5. Complete the final digitizing after receiving DSL responses;
6. Rerun `run_Tier3_4_final` or run `09_add_issue_IDs_to_set.ipynb` accordingly;
7. Update QAQC counts by excluding the completed records.

#### Report
The script `initial_status_report.ipynb` was initially applied to report the progress. 

## Notes

The scripts are organized to run by set. A set is collected when DSL sends data at certain time following their master schedule. The WD project "WD2019-0259" is from the non-participating counties. GIS files for the WD projects "WD2017-0229" and "WD2021-0179" are provided.

## WD SA data description

### Objectives and processes

DSL requested data attritutes and test data for their database development tests on June 2023. The process involves discussion on the sample data, steps to review study area polygons to get sample data, and follow-up on the criteria for final deliverables.

### Study area polygons

Run `join_WD_with_SA_by_taxmap` in `10_review_data_from_all_sets` to get the deliverable.