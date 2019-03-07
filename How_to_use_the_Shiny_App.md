
# Instruction to use the Shiny App:

Please read this document to use the [Shiny Web App for Data analysis tool and Sample Size Calculator for Non-Inferiority and Equivalence tests in SMART.](http://13.250.172.122/shiny/NI_EQ/)



## What can be done in the Shiny App:

1. Calculate the minimum sample size required for the Non-Inferiority test in SMART.

2. Calculate the minimum sample size required for the Equivalence test in SMART. 

3. Do the data analysis based on your completed SMART data.




## How to calculate the minimum sample size required for the Non-Inferiority test in SMART:

In a non-inferiority testing, the goal is to establish that a new intervention (say AI #1, which is cost-effective or less burdensome) yields favorable outcomes that, when compared to another active control (say AI #3) intervention (i.e., an intervention approach with established evidence of effectiveness), is not below some pre-stated non-inferiority margin.

1. Go to the [app.](http://13.250.172.122/shiny/NI_EQ/)
2. Select: Sample Size for Non Inferiority
3. Set the significance level. The default value is 0.05.
4. Set the desired power for the non-inferiority test. The default value is 0.80.
5. Set the standardized effect size. The default value is 0.3. The Standardized Effect Size = (The standardized Non-Inferiority margin) - (The standardized mean difference). See the 'Non-Inferiority Test' section of the paper given in reference for details.
6. The app will automatically update to show the minimum sample size for your inputed values in the last line of the page. For selected default values of above input parameters, the estimated total sample size (includes all AIs) needed in the SMART is 138.


## How to calculate the minimum sample size required for the Equivalence test in SMART:

In an equivalence testing, the goal is to establish that a new intervention (say AI #1, which is cost-effective or less burdensome) yields favorable outcomes that, when compared to another active control (say AI #3) intervention (i.e., an intervention approach with established evidence of effectiveness), is neither superior nor inferior with respect to a pre-specified equivalence margin.


1. Go to the [app.](http://13.250.172.122/shiny/NI_EQ/)
2. Select: Sample Size for Equivalence
3. Set the significance level. The default value is 0.05.
4. Set the standardized equivalence margin. The default value is 0.25.
5. Set the standardized effect size. The default value is 0.05. See the 'Equivalence Test' section of the paper given in reference for details.
6. Set a sample size (any approximate value) to start. The default value is 300. 
7. The app will automatically update to show the estimated power for selected sample size and the esitimated powers for next 20 values of the sample sample size. For default value, the app will show the powers corresponding to sample sizes 300-320.
8. Choose the sample size which shows the power that you are targetting. For example if you want to achieve 80% power then choose a sample size that gives atleast 0.80 power.
9. If there is no power close to 0.8 then repeat steps 6-8 with a different sample size (increase or decrease accrodingly).


## How to do the data analysis based on your completed SMART data:




## Reference:

Palash Ghosh, Inbal Nahum-Shani, Bonnie Spring and Bibhas Chakraborty (2018). Non-Inferiority and Equivalence Tests in A Sequential Multiple-Assignment Randomized Trial (SMART). Submitted.
