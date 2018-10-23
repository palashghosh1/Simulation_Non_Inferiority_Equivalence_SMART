

  #   title: "Simulation for Non-Inferiority test for SMART"
  #   author: "Dr. Palash Ghosh"
  #   date: "29/07/2016"
  #   last modified: 9 October 2018
  #   Description: This program calculates two differnet sample sizes (N) and powers based on 1000 simulations
  #                based on exact variance formula.



  
  # simulation size:

  sim.size = 1000 #1000 
  
  
  # specify the parameters for responders and non-responders
  
  # for responders:
  
  res.0  <- 0.02  ;   res.1.A  <-  0.8   ;   res.1.B  <- 0.7
  
  # for non-responders:
  
  nres.0 <- 0.03  ;   nres.1.A <-  0.5  ;   nres.1.B <- 0.4
  
  
  # Other parameters:
  
  nres.2.AD <- 1.3 
  
  nres.2.BC <- 0.8
  
  
  # specify the mean and variance of latent variable for treatment A
  
  mu.L.A <- 2
  sig.L <- 0.2  # note that it sigma not sigma^2
  


  ## main function:
  
  main_func <- function(t, sigma, theta, g1, g2, alpha, beta, for_delta ){
    

    # this parameter will control the vlaue of delta
 
    # for distinct-path: scenarios 1-25 are allocated for distinct-path simulations
    
    if(t <= 25){ 
      nres.2.AC <- for_delta
      nres.2.BD <- 0.8
      }
    
    # for shared-path: scenarios 26-50 are allocated for shared-path simulations
    
    if((t > 25) & (t <= 50)){  
      nres.2.AC <- 0.3  
      nres.2.BD <- for_delta
      }
    
    
    # calculate the cut-off based on latent variable for A
    
    eta <- mu.L.A + sig.L*qnorm(1-g1)
    
    
    # get the value of mu.L.B: mean of latent variable for treatment B
    # variance is same 
    
    mu.L.B <- eta - sig.L*qnorm(1-g2)
    mu.L.B <- mu.L.B[[1]]
    
    
    # truncated mean for non-responder: in leg start with A
    
    uuu <- (eta - mu.L.A)/sig.L
    mu.L.A_NR <- mu.L.A - sig.L*(dnorm(uuu)/pnorm(uuu))
    
    
    # truncated mean for non-responder: in leg start with B
    
    www <- (eta - mu.L.B)/sig.L
    mu.L.B_NR <- mu.L.B - sig.L*(dnorm(www)/pnorm(www))
    
    
    # to used in variance calculations
    
    mu.AA <- round((res.0 + res.1.A*mu.L.A), 2)
    mu.AC <- round((nres.0 + nres.1.A*mu.L.A + nres.2.AC*mu.L.A_NR), 2)
    mu.AD <- round((nres.0 + nres.1.A*mu.L.A + nres.2.AD*mu.L.A_NR), 2)
    
    
    mu.BB <- round((res.0 + res.1.B*mu.L.B), 2)
    mu.BC <- round((nres.0 + nres.1.B*mu.L.B + nres.2.BC*mu.L.B_NR), 2)
    mu.BD <- round((nres.0 + nres.1.B*mu.L.B + nres.2.BD*mu.L.B_NR), 2)
    
    
    # calculate the value of delta.dist:
    
    avg.AC <- (res.0 + res.1.A*mu.L.A)*g1 + (nres.0 + nres.1.A*mu.L.A + nres.2.AC*mu.L.A_NR)*(1-g1)
    avg.BC <- (res.0 + res.1.B*mu.L.B)*g2 + (nres.0 + nres.1.B*mu.L.B + nres.2.BC*mu.L.B_NR)*(1-g2)
    
    
    # calculate the effect size delta: to change the delta velue change the intail parameters:
    
    delta.dist <- round((avg.BC - avg.AC), 2)
    
    
    # calculate the value of delta.shared:

    avg.BC <- (res.0 + res.1.B*mu.L.B)*g2 + (nres.0 + nres.1.B*mu.L.B + nres.2.BC*mu.L.B_NR)*(1-g2)
    avg.BD <- (res.0 + res.1.B*mu.L.B)*g2 + (nres.0 + nres.1.B*mu.L.B + nres.2.BD*mu.L.B_NR)*(1-g2)
    
    # calculate the effect size delta: to change the delta velue change the intail parameters:

    delta.shared <- round((avg.BC - avg.BD), 2)
    
    

    # set z-alpha and z-beta of standard normal distribution:
    
    z_alpha <- qnorm(1 - alpha) ; z_beta <- qnorm(1 - beta)
    
    
    ## sample size formula:
    
    if(t == 1){print("Distinct-Path")}
    
    if(t <= 25){  
      
      # for distinct-path
     
      N.ext <- cal_N(sigma, theta, delta.dist, g1, g2, alpha, beta, distinct_path="Y", 
                 mu.AA, mu.AC, mu.AD, mu.BB, mu.BC, mu.BD, indi.noninf.equiv = "NI")
    }
    
    if(t == 26){print("Shared-Path")}
    
    if((t > 25) & (t <= 50)){  
     
      # for shared-path
      
      N.ext <- cal_N(sigma, theta, delta.shared, g1, g2, alpha, beta, distinct_path="N", 
                 mu.AA, mu.AC, mu.AD, mu.BB, mu.BC, mu.BD, indi.noninf.equiv = "NI")
    }
    
    if(N.ext > 10^5){ print("N is larger than 10^5. N has been set to 10^5")
      N.ext = 10^5
      print(N.ext)}
    
    print(c(t, N.ext))
    
    ## data generation for power calculation:
    
    n1.ext= round((N.ext/4)*(1+g1))
    n2.ext= round((N.ext/4)*(1+g2))
    
    # calculate variance for distinct-path
    
    V.sigma.DP.ext = ( 2*(4 - g1 -g2)*(sigma^2) + (  g1*(2-g1)*(mu.AA^2) + (3-2*g1-g1^2)*(mu.AC^2) -2*g1*(1-g1)*mu.AA*mu.AC
                                              + g2*(2-g2)*(mu.BB^2) + (3-2*g2-g2^2)*(mu.BC^2) -2*g2*(1-g2)*mu.BB*mu.BC ) )
    


    # calculate variance for shared-path
    
    V.sigma.SP.ext = ( 2*(4 - g2 -g2)*(sigma^2) + (  g2*(2-g2)*(mu.BB^2) + (3-2*g2-g2^2)*(mu.BC^2) -2*g2*(1-g2)*mu.BB*mu.BC
                                                   + g2*(2-g2)*(mu.BB^2) + (3-2*g2-g2^2)*(mu.BD^2) -2*g2*(1-g2)*mu.BB*mu.BD ) 
                       - 2*(2*g2*(sigma^2 + mu.BB^2) - (g2*mu.BB + (1-g2)*mu.BC)*(g2*mu.BB + (1-g2)*mu.BD)) )
    


    
    ## Power calculation:    
    # calling repeat function:
    
    # calling parallel::mclapply for program running on linux 
    
    if(Sys.info()[[1]] == "Linux"){
    
    II.ext <- unlist(parallel::mclapply(1:sim.size, function(i) func_repeat(i, N.ext, sigma, theta, g1, g2, alpha, beta, z_alpha, n1.ext, n2.ext,
                                                                        mu.L.A, sig.L,  res.0, res.1.A, res.1.B, nres.0, nres.1.A, nres.1.B, 
                                                                        nres.2.AD, nres.2.BC, nres.2.BD, nres.2.AC, 
                                                                        mu.AA, mu.AC, mu.AD, mu.BB, mu.BC, 
                                                                        seed = (1000+i)), mc.cores = detectCores() ))
    
    }
    
    
    if(Sys.info()[[1]] != "Linux"){
      
      II.ext <- unlist(lapply(1:sim.size, function(i) func_repeat(i, N.ext, sigma, theta, g1, g2, alpha, beta, z_alpha, n1.ext, n2.ext,
                                                                              mu.L.A, sig.L,  res.0, res.1.A, res.1.B, nres.0, nres.1.A, nres.1.B, 
                                                                              nres.2.AD, nres.2.BC, nres.2.BD, nres.2.AC, 
                                                                              mu.AA, mu.AC, mu.AD, mu.BB, mu.BC, 
                                                                              seed = (1000+i))))
      
    }
                                                                        
    
    II.ext <- matrix(II.ext, nrow=sim.size, byrow = T)
    
    
    # individual power:
    
    power.Di.NI.ext = mean(II.ext[,1])
    power.Sh.NI.ext = mean(II.ext[,4])
    power.Di.EQ.ext = mean(II.ext[,7])
    power.Sh.EQ.ext = mean(II.ext[,10])
    

    # creating output 
    output = matrix(0, 36, 1)
    
    output[1,1] = 0 
    output[2,1] = sigma
    output[3,1] = theta
    output[4,1] = delta.dist
    output[5,1] = delta.shared
    
    output[6,1] = alpha
    output[7,1] = 1 - beta
    output[8,1] = g1
    output[9,1] = g2

    # from exact variance:
    
    output[10,1] = N.ext
    output[11,1] = power.Di.NI.ext
    output[12,1] = power.Sh.NI.ext
    output[13,1] = power.Di.EQ.ext
    output[14,1] = power.Sh.EQ.ext
    
    # standadized effect size:
    
    output[15,1] = round((theta - delta.dist)/sqrt(V.sigma.DP.ext/2), 3)
    output[16,1] = round((delta.dist)/sqrt(V.sigma.DP.ext/2), 3)
    output[17,1] = round((theta - delta.shared)/sqrt(V.sigma.SP.ext/2), 3)   ## CHeck it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    output[18,1] = round((delta.shared)/sqrt(V.sigma.SP.ext/2), 3)
    

    return(t(output))
    
  }   # end of main function
  
  
  out <- matrix(0, 90, 36)


  ## calling the main function: t is the idex of the code-run
  
  ######## for 80% power:

  # for non-inferiority: calling the main_prog for distinct-path:
  

  out[1,] <- main_func(1, 3, 3,  0.3,  0.50, 0.05, 0.2, 0.01)
  out[2,] <- main_func(2, 3, 3,  0.3,  0.45, 0.05, 0.2, 0.01)
  out[3,] <- main_func(3, 3, 3,  0.3,  0.40, 0.05, 0.2, 0.01)
  out[4,] <- main_func(4, 3, 3,  0.3,  0.35, 0.05, 0.2, 0.01)
  out[5,] <- main_func(5, 3, 3,  0.3,  0.30, 0.05, 0.2, 0.01)

  out[6,] <- main_func(6, 3, 2.5,  0.3,  0.50, 0.05, 0.2, -0.2)
  out[7,] <- main_func(7, 3, 2.5,  0.3,  0.45, 0.05, 0.2, -0.2)
  out[8,] <- main_func(8, 3, 2.5,  0.3,  0.40, 0.05, 0.2, -0.2)
  out[9,] <- main_func(9, 3, 2.5,  0.3,  0.35, 0.05, 0.2, -0.2)
  out[10,] <- main_func(10, 3, 2.5,  0.3,  0.30, 0.05, 0.2, -0.2)
  
  
  
  

  # for non-inferiority: calling the main_prog for shared-path:


  out[26,] <- main_func(26, 3, 3,  0.30,  0.50, 0.05, 0.2, -0.38)
  out[27,] <- main_func(27, 3, 3,  0.30,  0.45, 0.05, 0.2, -0.38)
  out[28,] <- main_func(28, 3, 3,  0.30,  0.40, 0.05, 0.2, -0.38)
  out[29,] <- main_func(29, 3, 3,  0.30,  0.35, 0.05, 0.2, -0.38)
  out[30,] <- main_func(30, 3, 3,  0.30,  0.30, 0.05, 0.2, -0.38)

  out[31,] <- main_func(31, 3, 2.5,  0.30,  0.50, 0.05, 0.2, -0.53)
  out[32,] <- main_func(32, 3, 2.5,  0.30,  0.45, 0.05, 0.2, -0.53)
  out[33,] <- main_func(33, 3, 2.5,  0.30,  0.40, 0.05, 0.2, -0.53)
  out[34,] <- main_func(34, 3, 2.5,  0.30,  0.35, 0.05, 0.2, -0.53)
  out[35,] <- main_func(35, 3, 2.5,  0.30,  0.30, 0.05, 0.2, -0.53)
  
   
 
  # Creating output to print: 
  
  # for distinct-path
  
  out.print.Di <- round(out[1:10, c(2:3, 8:9, 10:11, 15)], 3)
  colnames(out.print.Di) <- c("sigma", "theta", 
                              "gamma.a", "gamma.ac",
                              "N", "Power.NI.Distinct-path", 
                              "Standard-Effect-Size")
                       
  

  # for shared-path
  
  out.print.Sh <- round(out[26:35, c(2:3, 8:9, 10, 12, 17)], 3)
  colnames(out.print.Sh) <- c("sigma", "theta", 
                              "gamma.a", "gamma.ac",
                              "N", 
                              "Power.NI.Shared-path", 
                              "Standard-Effect-Size")


  
  # exporting pdf files:
  
  
  pdf(file = "Distinct_path_NI_results.pdf", width = 30, height = 15)

  grid.table(out.print.Di)

  dev.off()
  
  pdf(file = "Shared_path_NI_results.pdf", width = 30, height = 15)
  
  grid.table(out.print.Sh)
  
  dev.off()


  
  
  
