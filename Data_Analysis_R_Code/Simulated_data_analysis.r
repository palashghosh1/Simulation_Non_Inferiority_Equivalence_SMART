
  #   title: "Data Analysis for Non-Inferiority and Equivalence tests for SMART"
  #   author: "Dr. Palash Ghosh"
  #   date: "29/10/2018"
  #   last modified: 26 November 2018
  #   Description: This program will read a SMART data set (with specific format and columns names)
  #                and then analyze the data set to compare two specific regimes in a non-inferiroity and/or 
  #                equivalence tests.
  
  
  # attaching the required libraries: 
  # If you do not have the library installed in your computer/laptop then please install them
  
  library(dplyr)
  library(readr)
  library(gridExtra)

  
  # attaching the different functions to be used
 
  source("test_power_data_analysis.r")


  
  
  
  # function to process the raw SMART data:


  function.data.analysis <- function(theta, alpha, treatment.seq, placebo.seq,
                                                        indi.shared.dist, indi.noninf.equiv){
    
    
    
    # reading the raw data:
    
    # Non-inferiority for distinct-path
    
    if (indi.shared.dist == "Di" & indi.noninf.equiv == "NI"){
    
    dat <- read.csv("simulated_SMART_data_Di_NI.csv")
    }
    
    
    # Non-inferiority for shared-path
    
    if (indi.shared.dist == "Sh" & indi.noninf.equiv == "NI"){
      
      dat <- read.csv("simulated_SMART_data_Sh_NI.csv")
    }
    
    
    
    
    # Equivalence for distinct-path
    
    if (indi.shared.dist == "Di" & indi.noninf.equiv == "EQ"){
      
      dat <- read.csv("simulated_SMART_data_Di_EQ.csv")
    }
    
    
    # Equivalence for shared-path
    
    if (indi.shared.dist == "Sh" & indi.noninf.equiv == "EQ"){
      
      dat <- read.csv("simulated_SMART_data_Sh_EQ.csv")
    }
    
    
    
    # making first and second stage treatments same for responders:
    
    # dat$Trt.2 <- ifelse(is.na(dat$Trt.2) & dat$R == 1, paste(dat$Trt.1), paste(dat$Trt.2))
    
    
    # make treatments and responder status factor:
    
    cols <- c("Trt.1", "R", "Trt.2")
    
    dat[cols] <- lapply(dat[cols], factor)
    
    
    
    # View(simulated_SMART_data)
    

    
    
  # Checking the Data:
    
    if (!(dim(dat)[2] == 5)) {
      
      cat("The data should contain exactly 5 columns", "\n") 
      
      } else {
        
        # cat("within 2nd loop", "\n")
        
        if(sum(colnames(dat) == c("sl", "Trt.1", "R", "Trt.2", "Y")) != 5){

          stop("The data should have exactly five columns with the specific order and names: sl, Trt.1, R, Trt.2, Y. Where 
            sl: serial numbers, Trt.1: First stage treatment, R: Response indicator, 
               Trt.2: Second stage treatment, Y: Continuous primary outcome")
          
        
        }}
    
    if(indi.shared.dist == "Di"){
      
      if((treatment.seq[1] == placebo.seq[1]))
      stop("First stage treatments in two distinct-path Adaptive Interventions should be different.")
    }
    
    if(indi.shared.dist == "Sh"){
      
      if((treatment.seq[1] != placebo.seq[1]))
        stop("First stage treatments in two shared-path Adaptive Interventions should be same.")
    }
    
      
      
    
    

    
    
    
  # alpha is the type-I error rate:
    
  z_alpha <- qnorm(1 - alpha)
  
  
  
  # N is the total number of observations:
  
  N <- nrow(dat)
  
  
  
  cat("The total number of observation in the SMART data:", N, "\n")
  cat("\n")
  

  
  
  # Observation per arm of the SMART:
  
  obs.per.arm <-  group_by(dat, Trt.1, Trt.2, R)%>%
                  summarise(N.obs = n())
  
  
  cat("Number of observations per arm of the SMART data: ", "\n")
  print(as.data.frame(obs.per.arm))
  
  
  cat("\n")
  cat("Treatment Sequence: ", "(", treatment.seq[1], ", ", treatment.seq[1], "^{R=1}", 
                                  treatment.seq[2], "^{R=0}", ")", sep="", "\n")
  cat("Placebo Sequence: ", "(", placebo.seq[1], ", ", placebo.seq[1], "^{R=1}", 
                              placebo.seq[2], "^{R=0}", ")", sep="", "\n")
  cat("\n")
  
  # cat(expression(paste("Cu"^"2+","at EC50",sep="")))
    
  
  
  # Observation per arm of the SMART:
  
  obs.1st.stage <- group_by(dat, Trt.1)%>%
                   summarise(N.Trt.1 = n())
  
  
  # Joined the above two data sets: obs.per.arm and obs.1st.stage
  
  joined.obs <- left_join(obs.per.arm, obs.1st.stage, by="Trt.1")
  
  
  
  # Calculate the response rates:
  
  res.rates.data <- filter(joined.obs, R == 1 )%>%
                    mutate(res.rate = round(N.obs/N.Trt.1,2))
  
  
  # Estimated g1 and g2 from data:
  
  g.treat <- res.rates.data$res.rate[res.rates.data$Trt.1 == treatment.seq[1]]
  g.placebo <- res.rates.data$res.rate[res.rates.data$Trt.1 == placebo.seq[1]]
  
  

  simu.out <- list()
  
  # calling the function: 
  # to_test <- function(data, treatment.seq, placebo.seq, indi.shared.dist, indi.noninf.equiv, N, g1.est, g2.est,
  # theta, z_alpha )
  
 
  power.out <- to_test(dat, treatment.seq, placebo.seq, indi.shared.dist, indi.noninf.equiv, 
                          N, g.treat, g.placebo, theta, z_alpha)
  

  
  # simu.out$I.Di.NI <- power.out$I
  # simu.out$Y.treat.avg.Di.NI <- power.out$Y.treat.avg
  # simu.out$Y.placebo.avg.Di.NI <- power.out$Y.placebo.avg
  
  simu.out$NI_EQ <- power.out$NI_EQ
  simu.out$distint.shared.path <- power.out$distint.shared.path
  
  simu.out$N <- N
  
  simu.out$Y.treat.avg.Di.NI <- power.out$Y.treat.avg
  simu.out$Y.placebo.avg.Di.NI <- power.out$Y.placebo.avg
 
  simu.out$p.val.Non_Inferiority <- power.out$p.val.Non_Inferiority
  simu.out$upp.bound.Bayes.fac.NI <- power.out$upp.bound.Bayes.fac.NI
    
  simu.out$p.val.Superiority <- power.out$p.val.Superiority
  simu.out$upp.bound.Bayes.fac.Sup <- power.out$upp.bound.Bayes.fac.Sup
  
  simu.out$Decision <- power.out$Decision
  
  
  

  return(simu.out)
  
  
  }
  
  
  # calling the function.data.analysis function:
  
  
  out.Di.NI <- function.data.analysis(theta = 2.5, alpha = 0.05, c("A", "C"), c("B", "C"),
                         indi.shared.dist = "Di", indi.noninf.equiv = "NI")
  print(t(as.data.frame(out.Di.NI)))
  cat("\n")
  cat("\n")
  
  

  out.Sh.NI <- function.data.analysis(theta = 2.5, alpha = 0.05, c("B", "D"), c("B", "C"),
                        indi.shared.dist = "Sh", indi.noninf.equiv = "NI")
  print(t(as.data.frame(out.Sh.NI)))
  
  cat("\n")
  cat("\n")
  


  out.Di.EQ <- function.data.analysis(theta = 2, alpha = 0.05, c("A", "C"), c("B", "C"),
                         indi.shared.dist = "Di", indi.noninf.equiv = "EQ")
  print(t(as.data.frame(out.Di.EQ)))
  cat("\n")
  cat("\n")
  
  
  
  out.Sh.EQ <- function.data.analysis(theta = 2, alpha = 0.05, c("B", "D"), c("B", "C"),
                                indi.shared.dist = "Sh", indi.noninf.equiv = "EQ")
  print(t(as.data.frame(out.Sh.EQ)))
  
  
  
  out <- rbind(t(as.data.frame(out.Di.NI)), "", t(as.data.frame(out.Sh.NI)), "", 
                t(as.data.frame(out.Di.EQ)), "", t(as.data.frame(out.Sh.EQ)))
   
   
  
  pdf("Results.pdf", height=13, width=8.5)

  grid.table(out)
 
  dev.off()
  
  
  
  
  