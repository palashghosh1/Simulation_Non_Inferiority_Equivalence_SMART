# Data Analysis for Non-Inferiority and Equivalence tests for SMART

This program will read a SMART data set (with specific format and columns names) and then analyze the data set to compare two specific regimes in a non-inferiroity and/or equivalence tests. If you are not familier with the statistical software R, then please use the shiny app for your data analysis from ![HERE]{https://palash.shinyapps.io/NI_EQ/}


## Data Description:

The four example simulated data sets are in csv format. All of them have precisely five variables (columns) with the specific variable names: sl, Trt.1, R, Trt.2, Y; where sl: denotes serial number (ex: 1,2,3,...), Trt.1: denotes first stage treatments (must be A or B), R: denotes responder status (must be 1 or 0, 1: responder, 0: non-responder), Trt.2: denotes second stage treatments (must be C or D if first stage treatment was A; or E or F if first stage treatment was B; or E=C and F=D), Y: denotes the continuous primary outcome. 

### List of r-files:

1. Simulated_data_analysis.r
2. test_power_data_analysis.r

### Simuated Example data sets:

1. simulated_SMART_data_Di_NI.csv
2. simulated_SMART_data_Sh_NI.csv
3. simulated_SMART_data_Di_EQ.csv
4. simulated_SMART_data_Sh_EQ.csv


### How to run the data analysis:
Intall the necessary libraries (dplyr, readr and gridExtra) in your computer. Put all the simulated data sets and two above r-codes in the same folder. Run the 'Simulated_data_analysis.r' in R to get the results given in the manuscript.


## Data analysis for new data:

Prepare your SMART data following the instruction given in the above Data Description. Modify the 'Simulated_data_analysis.r' accordingly with you new data file name and run the same.


