

  #   title: "Power Curves for Non-Inferiority and Equivalence test for SMART"
  #   author: "Dr. Palash Ghosh"
  #   date: "29/07/2016"
  #   last modified: 23 Oct 2018
  #   Description: This program generates power curves for Non-Inferiority and Equivalence test
  #                for SMART based on 1000 simulations.



  # simulation size:
  
  sim.size = 1000 

  
  
  # specify the parameters for responders and non-responders
  
  # for responders:
  
  res.0  <- 0.02  ;   res.1.A  <-  0.5   ;   res.1.B  <- 0.8
  
  # for non-responders:
  
  nres.0 <- 0.03  ;   nres.1.A <-  0.25  ;   nres.1.B <- 0.5
  
  # Other parameters:
  
  nres.2.AC <- 1.3

  nres.2.AD <- 1.3 ; 

  nres.2.BD <- 0.8

  
  # specify the mean and variance of latent variable for treatment A
  
  mu.L.A <- 2
  sig.L <- 0.2  # note that it sigma not sigma^2
  

 

  ## main function:
  
  main_func <- function(t, sigma, theta, g1, g2, alpha, beta, for_delta ){
    

    # this parameter will control the vlaue of delta

    nres.2.BC <- for_delta
    

    
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
    
    
    # mu.max <- round(max(mu.AA, mu.AC, mu.BB, mu.BC), 2)
    # mu.min <- round(min(mu.AA, mu.AC, mu.BB, mu.BC), 2)
    
    
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
    
    
   
  # setting different sample size for power curve generation:
    
    if(t <= 25){  N <- 100 }
    if((t > 25) & (t <= 50)){  N <- 200 }
    if((t > 50) & (t <= 75)){  N <- 300 }
    if((t > 75) & (t <= 100)){  N <- 500 }
    

    print(c(t, N))
    
    ## data generation for power calculation:
    
    n1= round((N/4)*(1+g1))
    n2= round((N/4)*(1+g2))


    # calculate variance:
    

    V.sigma.DP.ext = ( 2*(4 - g1 -g2)*(sigma^2) + (  g1*(2-g1)*(mu.AA^2) + (3-2*g1-g1^2)*(mu.AC^2) -2*g1*(1-g1)*mu.AA*mu.AC
                                                     + g2*(2-g2)*(mu.BB^2) + (3-2*g2-g2^2)*(mu.BC^2) -2*g2*(1-g2)*mu.BB*mu.BC ) )
    
    

    V.sigma.SP.ext = ( 2*(4 - g2 -g2)*(sigma^2) + (  g2*(2-g2)*(mu.BB^2) + (3-2*g2-g2^2)*(mu.BC^2) -2*g2*(1-g2)*mu.BB*mu.BC
                                                     + g2*(2-g2)*(mu.BB^2) + (3-2*g2-g2^2)*(mu.BD^2) -2*g2*(1-g2)*mu.BB*mu.BD ) 
                       - 2*(2*g2*(sigma^2 + mu.BB^2) - (g2*mu.BB + (1-g2)*mu.BC)*(g2*mu.BB + (1-g2)*mu.BD)) )
    
    
    ## Power calculation:    
    # calling repeat function:
    
    # calling parallel::mclapply for program running on linux 
    
    if(Sys.info()[[1]] == "Linux"){
      
      II <- unlist(parallel::mclapply(1:sim.size, function(i) func_repeat(i, N, sigma, theta, g1, g2, alpha, beta, z_alpha, n1, n2,
                                                                              mu.L.A, sig.L,  res.0, res.1.A, res.1.B, nres.0, nres.1.A, nres.1.B, 
                                                                              nres.2.AD, nres.2.BC, nres.2.BD, nres.2.AC, 
                                                                              mu.AA, mu.AC, mu.AD, mu.BB, mu.BC, 
                                                                              seed = (100+i)), mc.cores = detectCores() ))
      
    }
    
    
    if(Sys.info()[[1]] != "Linux"){
      
      II <- unlist(lapply(1:sim.size, function(i) func_repeat(i, N, sigma, theta, g1, g2, alpha, beta, z_alpha, n1, n2,
                                                                  mu.L.A, sig.L,  res.0, res.1.A, res.1.B, nres.0, nres.1.A, nres.1.B, 
                                                                  nres.2.AD, nres.2.BC, nres.2.BD, nres.2.AC, 
                                                                  mu.AA, mu.AC, mu.AD, mu.BB, mu.BC, 
                                                                  seed = (100+i))))
      
    }
    

    II <- matrix(II, nrow=sim.size, byrow = T)
    
    
    # individual power:
    
    power.Di.NI = mean(II[,1])
    power.Sh.NI = mean(II[,4])
    power.Di.EQ = mean(II[,7])
    power.Sh.EQ = mean(II[,10])


    
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
    output[10,1] = n1 
    output[11,1] = n2 
    output[12,1] = N

    output[13,1] = power.Di.NI
    output[14,1] = power.Sh.NI
    output[15,1] = power.Di.EQ
    output[16,1] = power.Sh.EQ

    output[17,1] = avg.AC
    output[18,1] = avg.BC
    output[19,1] = avg.BD

    output[20,1] = mean(II[,2])
    output[21,1] = mean(II[,3])
    # output[22,1] = mean(II[,6])
    output[22,1] = mean(II[,5])

    output[23,1] = nres.2.BC


    output[33,1] = round((theta - delta.dist)/sqrt(V.sigma.DP.ext/2), 5)
    output[34,1] = round((delta.dist)/sqrt(V.sigma.DP.ext/2), 5)
    output[35,1] = round((theta - delta.shared)/sqrt(V.sigma.SP.ext/2), 5)   
    output[36,1] = round((delta.shared)/sqrt(V.sigma.SP.ext/2), 5)

    
    return(t(output))
    
  }   # end of main function
  
  
  out <- matrix(0, 100, 36)


  ## calling the main function: t is the idex of the code-run
  

  out[1,] <- main_func(1, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.8)
  out[2,] <- main_func(2, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.5)
  out[3,] <- main_func(3, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.3)
  out[4,] <- main_func(4, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.0)
  out[5,] <- main_func(5, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.8)
  out[6,] <- main_func(6, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.6)
  out[7,] <- main_func(7, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.4)
  out[8,] <- main_func(8, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.2)
  out[9,] <- main_func(9, 2, 2,  0.3,  0.4, 0.05, 0.2, 0)
  out[10,] <- main_func(10, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.2)
  out[11,] <- main_func(11, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.4)
  out[12,] <- main_func(12, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.6)
  out[13,] <- main_func(13, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.8)
  out[14,] <- main_func(14, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.0)
  out[15,] <- main_func(15, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.2)
  out[16,] <- main_func(16, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.4)
  out[17,] <- main_func(17, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.6)
  out[18,] <- main_func(18, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.8)
  out[19,] <- main_func(19, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.0)
  out[20,] <- main_func(20, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.2)
  out[21,] <- main_func(21, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.4)
  out[22,] <- main_func(22, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.0)
  out[23,] <- main_func(23, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.2)
  out[24,] <- main_func(24, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.4)
  out[25,] <- main_func(25, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.6)
  
  
  

  out[26,] <- main_func(26, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.8)
  out[27,] <- main_func(27, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.5)
  out[28,] <- main_func(28, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.3)
  out[29,] <- main_func(29, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.0)
  out[30,] <- main_func(30, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.8)
  out[31,] <- main_func(31, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.6)
  out[32,] <- main_func(32, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.4)
  out[33,] <- main_func(33, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.2)
  out[34,] <- main_func(34, 2, 2,  0.3,  0.4, 0.05, 0.2, 0)
  out[35,] <- main_func(35, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.2)
  out[36,] <- main_func(36, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.4)
  out[37,] <- main_func(37, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.6)
  out[38,] <- main_func(38, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.8)
  out[39,] <- main_func(39, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.0)
  out[40,] <- main_func(40, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.2)
  out[41,] <- main_func(41, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.4)
  out[42,] <- main_func(42, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.6)
  out[43,] <- main_func(43, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.8)
  out[44,] <- main_func(44, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.0)
  out[45,] <- main_func(45, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.2)
  out[46,] <- main_func(46, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.4)
  out[47,] <- main_func(47, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.0)
  out[48,] <- main_func(48, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.2)
  out[49,] <- main_func(49, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.4)
  out[50,] <- main_func(50, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.6)


  out[51,] <- main_func(51, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.8)
  out[52,] <- main_func(52, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.5)
  out[53,] <- main_func(53, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.3)
  out[54,] <- main_func(54, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.0)
  out[55,] <- main_func(55, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.8)
  out[56,] <- main_func(56, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.6)
  out[57,] <- main_func(57, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.4)
  out[58,] <- main_func(58, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.2)
  out[59,] <- main_func(59, 2, 2,  0.3,  0.4, 0.05, 0.2, 0)
  out[60,] <- main_func(60, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.2)
  out[61,] <- main_func(61, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.4)
  out[62,] <- main_func(62, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.6)
  out[63,] <- main_func(63, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.8)
  out[64,] <- main_func(64, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.0)
  out[65,] <- main_func(65, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.2)
  out[66,] <- main_func(66, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.4)
  out[67,] <- main_func(67, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.6)
  out[68,] <- main_func(68, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.8)
  out[69,] <- main_func(69, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.0)
  out[70,] <- main_func(70, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.2)
  out[71,] <- main_func(71, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.4)
  out[72,] <- main_func(72, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.0)
  out[73,] <- main_func(73, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.2)
  out[74,] <- main_func(74, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.4)
  out[75,] <- main_func(75, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.6)

  out[76,] <- main_func(76, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.8)
  out[77,] <- main_func(77, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.5)
  out[78,] <- main_func(78, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.3)
  out[79,] <- main_func(79, 2, 2,  0.3,  0.4, 0.05, 0.2, 1.0)
  out[80,] <- main_func(80, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.8)
  out[81,] <- main_func(81, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.6)
  out[82,] <- main_func(82, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.4)
  out[83,] <- main_func(83, 2, 2,  0.3,  0.4, 0.05, 0.2, 0.2)
  out[84,] <- main_func(84, 2, 2,  0.3,  0.4, 0.05, 0.2, 0)
  out[85,] <- main_func(85, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.2)
  out[86,] <- main_func(86, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.4)
  out[87,] <- main_func(87, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.6)
  out[88,] <- main_func(88, 2, 2,  0.3,  0.4, 0.05, 0.2, -0.8)
  out[89,] <- main_func(89, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.0)
  out[90,] <- main_func(90, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.2)
  out[91,] <- main_func(91, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.4)
  out[92,] <- main_func(92, 2, 2,  0.3,  0.4, 0.05, 0.2, -2.6)
  out[93,] <- main_func(93, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.8)
  out[94,] <- main_func(94, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.0)
  out[95,] <- main_func(95, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.2)
  out[96,] <- main_func(96, 2, 2,  0.3,  0.4, 0.05, 0.2, 3.4)
  out[97,] <- main_func(97, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.0)
  out[98,] <- main_func(98, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.2)
  out[99,] <- main_func(99, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.4)
  out[100,] <- main_func(100, 2, 2,  0.3,  0.4, 0.05, 0.2, 2.6)
  

  
  
  
  
  dd <- out[1:100, ]

  d.Di.NI <- data.frame(effect_size = dd[ , 33], power = dd[ , 13], N = as.factor(dd[ , 12]))
  d.Sh.NI <- data.frame(effect_size = dd[ , 35], power = dd[ , 14], N = as.factor(dd[ , 12]))
  d.Di.EQ <- data.frame(effect_size = dd[ , 33], power = dd[ , 15], N = as.factor(dd[ , 12]))
  d.Sh.EQ <- data.frame(effect_size = dd[ , 35], power = dd[ , 16], N = as.factor(dd[ , 12]))
  
  
  # plots:
  
  
  plot.1 <- ggplot(data=d.Di.NI, aes(x=effect_size, y=power, group = N, color = N, linetype=N)) +
    geom_line( ) +
    ggtitle("Non-Inferiority: Distinct-path \n") + 
    xlab("Standardized effect size \n")  +
    ylab("Power \n") 
  
  
  plot.2 <- ggplot(data=d.Sh.NI, aes(x=effect_size, y=power, group = N, color = N, linetype=N)) +
    geom_line( ) +
    ggtitle("Non-Inferiority: Shared-path \n") + 
    xlab("Standardized effect size \n")  +
    ylab("Power \n") 
  
  
  plot.3 <- ggplot(data=d.Di.EQ, aes(x=effect_size, y=power, group = N, color = N, linetype=N)) +
    geom_line( ) +
    ggtitle("Equivalence: Distinct-path \n") + 
    xlab("Standardized effect size \n")  +
    ylab("Power \n") 
  
  
  plot.4 <- ggplot(data=d.Sh.EQ, aes(x=effect_size, y=power, group = N, color = N, linetype=N)) +
    geom_line( ) +
    ggtitle("Equivalence: Shared-path \n") + 
    xlab("Standardized effect size \n")  +
    ylab("Power \n") 
  
  multiplot(plot.1, plot.2, plot.3, plot.4, cols=2)
  
  
  pdf("Power_curve.pdf", width = 12, height = 8)
  multiplot(plot.1, plot.2, plot.3, plot.4, cols=2) # Make plot
  dev.off()
  
  # postscript("Power_curve.eps", width = 480, height = 480)
  # multiplot(plot.1, plot.2, plot.3, plot.4, cols=2) # Make plot
  # dev.off()
  
  jpeg("Power_curve.jpeg", width = 1024, height = 768)
  multiplot(plot.1, plot.2, plot.3, plot.4, cols=2) # Make plot
  dev.off()
  
  tiff("Power_curve.tiff", width = 1024, height = 768)
  multiplot(plot.1, plot.2, plot.3, plot.4, cols=2) # Make plot
  dev.off()
  
 
  

  
  
  