This repository contains scripts to recreate the analysis and figures in "Linking community-climate disequilibrium with ecosystem function" by Stemkovski et al. All data are generated in the scripts, so no data files are included.

Run the scripts in the following order:
1. main_figures.R
2. acclimation_model_tests.R
3. disequilibrium_response_forms.R
4. simulation_demo.R
5. apply_to_simulation.R

The following auxiliary scripts are called from those above: functions.R, CommunityTempDis.R, and SpeciesPoolGen.R. The latter 2 are found in the lib folder and are associated with the metacommunity simulation from Alexander et al. (2018) and Adler et al. (2020).

The following R-packages are required: data.table, fields.

Please email Michael Stemkovski at m.stemkovski@gmail.com with any questions.


