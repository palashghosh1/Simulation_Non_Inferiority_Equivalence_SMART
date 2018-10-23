

  # This function will call the non-inferiority, equivalence and power-curve main prog
  # This genereate all the outputs given in the article in pdf format
  # This will also generates power curve graph in pdf, tiff and jpeg formats


  # attaching the required libraries: 
  # If you do not have the library installed in your computer/laptop then please install them
  
  library(dplyr)
  library(ggplot2)
  library(xtable)
  library(parallel)
  library(gridExtra)
  
  
  # attaching the different functions to be used
  
  source("gen_data.r")  
  source("test_power.r")
  source("multiplot.r")
  
  
  # calling the non-inferiority simulation function:
  
  print("Simulations for Non-Inferiority test:")
  
  source("main_prog_NON_INFERIORITY.r")
  
  
  # calling the equivalence simulation function:
  
  print("Simulations for Equivalence test:")
  
  
  source("main_prog_EQUIVALENCE.r")
  
  
  # calling the power curve generation function:
  
  print("Power Curves:")
  
  
  source("main_prog_POWER_CURVE.r")
  
  