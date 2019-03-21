
# Instruction to use the Shiny App:

Please read this document to use the [Shiny Web App for Data analysis tool and Sample Size Calculator for Non-Inferiority and Equivalence tests in a SMART.](https://palash.shinyapps.io/NI_EQ/)



## What can you do using the Shiny App?:

1. Calculate the minimum sample size required for a Non-Inferiority test in a SMART.

2. Calculate the minimum sample size required for an Equivalence test in a SMART. 

3. Analyze data from your completed SMART.




## How to calculate the minimum sample size required for a Non-Inferiority test in a SMART?:

In a non-inferiority testing, the goal is to establish that a new adaptive intervention (say AI #1, which may be cost-effective or less burdensome) yields favorable outcomes that, when compared to an active control (say AI #3) intervention (i.e., an intervention approach with established evidence of effectiveness), is not below some pre-specified non-inferiority margin.

1. Go to the [app.](https://palash.shinyapps.io/NI_EQ/)
2. Select: Sample Size for Non-Inferiority
3. Set the significance level. The default value is 0.05.
4. Set the desired power for the non-inferiority test. The default value is 0.80.
5. Set the standardized non-inferiority margin. The default value is 0.5.
6. Set the standardized mean difference. The default value is 0.2. 
8. The Standardized Effect Size = (The standardized Non-Inferiority margin) - (The standardized mean difference). See the 'Non-Inferiority Test' section of the paper given in reference for details.
8. The webpage will automatically refresh to show the minimum sample size based on the inputted values in the last line of the page. For selected default values of the above input parameters, the estimated total sample size (includes all AIs) needed in the SMART is 138.


## How to calculate the minimum sample size required for an Equivalence test in a SMART?:

In an equivalence testing, the goal is to establish that a new intervention (say AI #1, which may be cost-effective or less burdensome) yields favorable outcomes that, when compared to an active control (say AI #3) intervention (i.e., an intervention approach with established evidence of effectiveness), is neither superior nor inferior with respect to a pre-specified equivalence margin.


1. Go to the [app.](https://palash.shinyapps.io/NI_EQ/)
2. Select: Sample Size for Equivalence
3. Set the significance level. The default value is 0.05.
4. Set the standardized equivalence margin. The default value is 0.25.
5. Set the standardized mean difference. The default value is 0.05. See the 'Equivalence Test' section of the paper given in reference for details.
6. Set a sample size (any approximate value) to start. The default value is 300. 
7. The webpage will automatically refresh to show the estimated power for the selected sample size and the estimated powers for next 20 values of the sample size. For the default value, the app will show the powers corresponding to sample sizes 300-320.
8. Choose the sample size which shows the power that you are targetting. For example, if you want to achieve 80% power then choose a sample size that gives at least 0.80 power.
9. If there is no power close to 0.80, then repeat steps 6-8 with a different sample size (increase or decrease accordingly).


## How to do the data analysis based on your completed SMART data?:


![shematic_diagram_ppt](https://user-images.githubusercontent.com/43629013/54010544-273bff80-41aa-11e9-880f-4970b0de6653.jpg)


The data analysis is based on the SMART depicted in the above schematic diagram. 

1. Go to the [app.](https://palash.shinyapps.io/NI_EQ/)
2. Select: Data Analysis
3. Prepare your input data: The input data must be in excel (xls or xlsx) format. It should contain precisely five variables (columns) with the specific variable names: sl, Trt.1, R, Trt.2, Y; where sl: denotes serial number (ex: 1,2,3,...), Trt.1: denotes first stage treatments (must be A or B), R: denotes responder status (must be 1 or 0; 1: responder, 0: non-responder), Trt.2: denotes second stage treatments (must be A, C or D if first stage treatment was A; or B, E or F if first stage treatment was B), Y: denotes the continuous primary outcome. 
4. Upload your data by clicking the "Browse" under the heading "Reading user's Data:". 
5. Set the significance level. The default value is 0.05.
6. Set the Non-Inferiority/Equivalence Margin. The default value is 2.5.
7. Choose the test for Data Analysis: Non Inferiority or Equivalence.
8. Choose the comparison type: Distinct-Path or Shared-Path.
9. Choose the first stage treatment for New AI: A or B.
10. Choose the second stage treatment for non-responders in New AI: C or D (if 'A' in step 9); E or F (if 'B' in step 9).
11. The first stage treatment for Control AI: A or B (Automatically chosen based on your previous choices).
12. Choose the Second Stage Treatment for non-responders in Control AI: from (E or F) or from (C or D) based on previous choices.
13. See the end section "Data Analysis Results" for detailed output. The result includes the p-values and the decision to reject or not to reject the null hypothesis. 



## Reference:

Palash Ghosh, Inbal Nahum-Shani, Bonnie Spring and Bibhas Chakraborty (2019). Non-Inferiority and Equivalence Tests in A Sequential Multiple-Assignment Randomized Trial (SMART). Under Revision.
