Scripts needed to make variant maps and overlap maps

-CreateVariantMaps.m takes spatial correlation maps and outputs variants with unique IDs
Spatial correlation maps are located in ~/Box/Quest_Backup/HCP_analyses/from_Ben. HCP_variants/spCorr. Currently uses two thresholds: 5% and 10%. Existing variant maps can be found in ~/Box/DATA/HCP/variantMaps.

-overlapmap.m takes variant maps with unique IDs, turns all unique IDs to 1's, and sums across subjects to determine the proportion of subjects that have a variant at any given vertex. Saves proportions to cifti.

-goodSubs_allInfo.xlsx provides the subjects IDs for the three groups: left handers, right handers, and 'middle' handers.

-makingdatafiles.m uses subject IDs found in goodSubs_allInfo.xlsx to create a variable with the paths to the subjects' files. Both scripts above use this function.

