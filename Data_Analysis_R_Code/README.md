# Data Analysis for Non-Inferiority and Equivalence tests for SMART

This program will read a SMART data set (with specific format and columns names) and then analyze the data set to compare two specific regimes in a non-inferiroity and/or equivalence tests.


## Data Description:

The input data must be csv format. It should contain precisely five variables (columns) with the specific variable names: sl, Trt.1, R, Trt.2, Y; where sl: denotes serial number (ex: 1,2,3,...), Trt.1: denotes first stage treatments (must be A or B), R: denotes responder status (must be 1 or 0, 1: responder, 0: non-responder), Trt.2: denotes second stage treatments (must be C or D if first stage treatment was A; or E or F if first stage treatment was B), Y: denotes the continuous primary outcome. 
