

# Data treatment files

This folder ("Data Treatment") is made for people using the PTR-MS for H2S measurements on AU. 
It makes data treatment easier, by using the raw excel file obtained from the PTR-MS directly, 
and exporting the time and concentration as a csv file.

# How to use the scripts 
To use the script, start på putting all your raw data (excel files) in the "Raw_data" folder. 
Also put your calibration file in that folder. 

## Calibration
A calibration file should be made by measuring a known concentration at a variety of water content / air humidity. 
This can eg. be done in a Tedlar bag that is injectred with water. 

Then, go to /Data treatment/Calibration_func, and adjust the known concentration to the one you measured at, adjust the start and 
stop cycles (rows in the excelfile) if you dont want all the data to be a part of the calibration. 
Also adjust the file name to whatever you called your calibration file. 
If you are in eg. Spyder, a graph with your calibration data and the obtained calibration curve will be shown. 
It is a good idea to check that this looks reasonable.

## Getting your data
After the calibration, go to /Data Tretmeant/DT_func" and change the filename to your first file. Run the file. Change the file name to
the next file. Repeat for all your files. The csv files should turn up in the "Processed_data" folder with the same 
names that they originally had. 

## Plotting your data 
I you are interested in this specific projects plotting, go to /scripts/Plotting_standard_experiments and adjust the number of experiment 
- note that this has to be done 2 places = 4 numbers need to be changed.This will overwrite the plots currently in the "Plots" folder.


# Creator and contact
This is still not tested by others than me, Anne Nørgaard Mortensen. Email: 201806874@post.au.dk
