
# Shiny APP:

[Shiny Web App for Sample Size Calculation and Data Analysis](http://13.250.172.122/shiny/NI_EQ/) 


# Simulations for Non-Inferiority and Equivalence Tests in A Sequential Multiple-Assignment Randomized Trial (SMART)

This project will generate all the simulation tables and power curves used in the manuscript
"Non-Inferiority and Equivalence Tests in A Sequential Multiple-Assignment Randomized Trial (SMART)" 
authored by Palash Ghosh, Inbal Nahum-Shani, Bonnie Spring and Bibhas Chakraborty.

## List of r-files:

Copy all the codes in your computer/laptop. The different files are:

1. call_main_progs.r
2. main_prog_NON_INFERIORITY.r
3. main_prog_EQUIVALENCE.r
4. main_prog_POWER_CURVE.r
5. gen_data.r
6. test_power.r
7. multiplot.r



### Prerequisites R-libraries

Install following R-libraries inside R-software in your computer/laptop:

  1. dplyr
  2. ggplot2
  3. xtable
  4. parallel
  5. gridExtra


### How to run the simulation codes

Once you have installed all the libraries in your computer/laptop, run the `call_main_progs.r' in your Rstudio/R. 
On completion, it will generate all the simulation tables and power curves used in the manuscript. 


## Operating system

The codes have been tested in both linux [Processor	: 32x Intel(R) Xeon(R) CPU E5-2665 0 @ 2.40GHz, 
Memory: 65890MB, Operating System	: Ubuntu 16.04.3 LTS] and windows (Windows 10) system. 

Note that, the parallel programming will only work in the Linux system.


## Authors

* **Dr. Palash Ghosh** 


