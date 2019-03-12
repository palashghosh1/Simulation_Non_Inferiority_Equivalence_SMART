# Data Analysis for Non-Inferiority and Equivalence tests for SMART

This R-program will read a SMART data set (with specific format and columns names) and then analyze the data set to compare two specific regimes in a non-inferiroity and/or equivalence tests. If you are not familier with the statistical software R, then please use the shiny app for your data analysis from [HERE](https://palash.shinyapps.io/NI_EQ/). The insitructions to use the shiny app are given [HERE](https://github.com/palashghosh1/Simulation_Non_Inferiority_Equivalence_SMART/blob/master/How_to_use_the_Shiny_App.md)


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

1. Prepare your SMART data following the instructions given in the above Data Description.
2. Put the above two r-files (Simulated_data_analysis.r and test_power_data_analysis.r) and your new SMART data set in a same folder.
3. Open the 'Simulated_data_analysis.r'.
4. Inside the function 'function.data.analysis' in the above r-file, replace all the four simulated example datasets (in ".csv" format) by your new SMART data set name (in csv format).
5. Go to the end of the r-file at 'calling the function.data.analysis function:'. 
6. There are four examples given to call the 'function.data.analysis': 
7. First one, "out.Di.NI": Non-Inferiority test to compare two distinct-path AIs. 
8. Second one, "out.Sh.NI": Non-Inferiority test to compare two shared-path AIs. 
9. Third one, "out.Di.EQ": Equivalence test to compare two distinct-path AIs. 
10. Fourth one. "out.Sh.EQ": Equivalence test to compare two shared-path AIs. 


