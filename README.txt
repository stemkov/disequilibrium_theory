This repository contains scripts to recreate the analysis and figures in "Linking community-climate disequilibrium with ecosystem function" by Stemkovski et al. All data are generated in the scripts, so no data files are included.

Run the scripts in the following order:
1. main_figures.R
2. acclimation_model_tests.R
3. disequilibrium_response_forms.R
4. simulation_demo.R

The following auxiliary scripts are called from those above: functions.R, CommunityTempDis.R, and SpeciesPoolGen.R. The latter 2 are found in the "lib" folder and are associated with the metacommunity simulation from Alexander et al. (2018) and Adler et al. (2020).

Matlab scripts for reproducing figures from the supplement at in the "sup" folder. See README_sup.txt in the "sup" folder for more information about the scripts.

The following R-packages are required: data.table, fields.

You can also interact with a version of the figures from main_figures.R in an accompanying R-Shiny app here: https://www.stemkovski.com/visualizations/ecological-acclimation-app

Please email Michael Stemkovski at m.stemkovski@gmail.com with any questions.


