#!/bin/bash
myArray=(1a_fatdepot_import.R 2a_fatdepot_analyse.R 3a_fatdepot_additionalanalysis.R 4a_fatdepot_plot.R 5a_fatdepot_writetables.R 6a_fatdepot_writetext.R)
for str in ${myArray[@]}; do
chmod u+x ./$str
done
echo 'Initializing 1a_fatdepot_import.R' && ./1a_fatdepot_import.R && echo 'Initializing 2a_fatdepot_analyse.R' && ./2a_fatdepot_analyse.R && echo 'Initializing 3a_fatdepot_additionalanalysis.R' && ./3a_fatdepot_additionalanalysis.R && echo 'Initializing 4a_fatdepot_plot.R' && ./4a_fatdepot_plot.R && echo 'Initializing 5a_fatdepot_writetables.R' && ./5a_fatdepot_writetables.R && echo 'Initializing 6a_fatdepot_writetext.R' && ./6a_fatdepot_writetext.R && echo 'The master script finished without errors'
