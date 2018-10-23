

  # First function will generate sample size and
  # the second function will generate the data for non-inferiority and equivalence SMART:
  
  
  # calculate the sample size:

  cal_N <- function(sigma, theta, delta, g1, g2, alpha, beta, distinct_path, mu.AA, mu.AC, mu.AD, mu.BB, mu.BC, mu.BD, indi.noninf.equiv){
    
    z_alpha = qnorm(1 - alpha) ; z_beta = qnorm(1 - beta)
    z_beta.half = qnorm(1 - (beta/2))  # required for equivalence test
    
    # sample size formula:
    if(distinct_path == "Y"){
      
      V.sigma = ( 2*(4 - g1 -g2)*(sigma^2) + (  g1*(2-g1)*(mu.AA^2) + (3-2*g1-g1^2)*(mu.AC^2) -2*g1*(1-g1)*mu.AA*mu.AC
                                                      + g2*(2-g2)*(mu.BB^2) + (3-2*g2-g2^2)*(mu.BC^2) -2*g2*(1-g2)*mu.BB*mu.BC ) )

      
      

      if(indi.noninf.equiv == "NI"){
        
      std.eff.sz <- round( ((theta - delta)/sqrt(V.sigma/2)), 3)
      
      N = ceiling(2*((z_alpha + z_beta)^2)*(1/(std.eff.sz^2)) )
      
      }
      
      if(indi.noninf.equiv == "EQ"){
        
      std.eff.sz <- round( ((theta)/sqrt(V.sigma/2)), 3)

      N = ceiling(2*((z_alpha + z_beta.half)^2)*(1/(std.eff.sz^2)) )
        
      }

    }

    if(distinct_path == "N"){
      
        V.sigma = ( 2*(4 - g2 -g2)*(sigma^2) + (  g2*(2-g2)*(mu.BB^2) + (3-2*g2-g2^2)*(mu.BC^2) -2*g2*(1-g2)*mu.BB*mu.BC
                                                + g2*(2-g2)*(mu.BB^2) + (3-2*g2-g2^2)*(mu.BD^2) -2*g2*(1-g2)*mu.BB*mu.BD ) 
                    - 2*(2*g2*(sigma^2 + mu.BB^2) - (g2*mu.BB + (1-g2)*mu.BC)*(g2*mu.BB + (1-g2)*mu.BD)) )

      

      if(indi.noninf.equiv == "NI"){
        
        std.eff.sz <- round( ((theta - delta)/sqrt(V.sigma/2)), 3)
        
        N = ceiling(2*((z_alpha + z_beta)^2)*(1/(std.eff.sz^2)) )
        
      }
      
      if(indi.noninf.equiv == "EQ"){
        
        std.eff.sz <- round( ((theta)/sqrt(V.sigma/2)), 3)
        
        N = ceiling(2*((z_alpha + z_beta.half)^2)*(1/(std.eff.sz^2)) )
        
      }

      }
    
    if(!is.numeric(N)){ print("Unable to calculate N") }
    
    return(N)
  }


  # function to do the simulation  
  
  func_repeat <- function(i, N, sigma, theta, g1, g2, alpha, beta, z_alpha, n1, n2,
                        mu.L.A, sig.L,  res.0, res.1.A, res.1.B, nres.0, nres.1.A, nres.1.B, 
                        nres.2.AD, nres.2.BC, nres.2.BD, nres.2.AC, 
                        mu.AA, mu.AC, mu.AD, mu.BB, mu.BC, 
                        seed){
  
  set.seed(seed)
  

  # generate latent variable for treatment A: half of N patients get treatment A
  
  N.A <- ceiling(N/2)
  N.B <- N - N.A
  L.A <- round(rnorm(N.A, mu.L.A, sig.L), 2)
  
  
  # now estimate all the population parameter based on estimated g1 and g2:
  
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
  

  # generate the latent variable for treatment B:
  
  L.B <- round(rnorm(N.B, mu.L.B, sig.L), 2)
  

  
  # construct the data:
  
  dat <- data.frame(sl = 1:N, Latent.1 = c(L.A, L.B)) 
  
  dat <- mutate(dat, Trt.1 = ifelse(sl <= N.A, "A", "B"),        # allocate the 1st stage treatments
                R  = ifelse(Latent.1 > eta, 1, 0))%>%      # define the responder and non-responder based on cut-off value eta
    group_by(Trt.1, R)%>%
    mutate(Lat.1.mu = round(mean(Latent.1), 2),
                           # I = rbinom(n(), 1, 0.5),
           I = sample(c(rep(0, ceiling(n()/2)), rep(1, n() - ceiling(n()/2))) ),
           Trt.2 = ifelse(I == 1, "C", "D"),
           Trt.2 = ifelse(R == 0,  Trt.2, NA_character_))
  

  
  # generate mean and SD for primary outcome generation:
  
  
  dat <- group_by(dat, Trt.1, R, Trt.2)%>%
    mutate( Latent.2 = as.numeric(ifelse( (R==0 & Trt.1=="A" & Trt.2=="C"), round(rnorm(n(), mu.L.A_NR, sig.L), 2),
                                          ifelse( (R==0 & Trt.1=="A" & Trt.2=="D"), round(rnorm(n(), mu.L.A_NR, sig.L), 2),
                                                  ifelse( (R==0 & Trt.1=="B" & Trt.2=="C"), round(rnorm(n(), mu.L.B_NR, sig.L), 2),
                                                          ifelse( (R==0 & Trt.1=="B" & Trt.2=="D"), round(rnorm(n(), mu.L.B_NR, sig.L), 2), NA_real_))))),
            
            Lat.2.mu = ifelse(R == 0, round(mean(Latent.2), 2), NA_real_),
            
            mu.y = ifelse( (R == 1 & Trt.1=="A"), res.0 + res.1.A*Lat.1.mu, 
                           ifelse( (R == 1 & Trt.1=="B"), res.0 + res.1.B*Lat.1.mu,
                                   ifelse( (R == 0 & Trt.1=="A" & Trt.2=="C"),  nres.0 + nres.1.A*Lat.1.mu + nres.2.AC*Lat.2.mu,
                                           ifelse( (R == 0 & Trt.1=="A" & Trt.2=="D"),  nres.0 + nres.1.A*Lat.1.mu + nres.2.AD*Lat.2.mu,
                                                   ifelse( (R == 0 & Trt.1=="B" & Trt.2=="C"), nres.0 + nres.1.B*Lat.1.mu + nres.2.BC*Lat.2.mu,
                                                           nres.0 + nres.1.B*Lat.1.mu + nres.2.BD*Lat.2.mu )))) ) )
  
  
  
  
  
  data.leg <- group_by(dat, Trt.1, R, Trt.2)%>%
    summarise( mu.leg = mean(mu.y),
               p.leg = round(n()/N, 5),
               regime.n = n(),
               Lat.2.mu = mean(Latent.2))
  
  
  mu.leg <- data.leg$mu.leg
  p.leg  <- data.leg$p.leg
 
  sigma.star  <- sigma
  
  dat <- dat%>%group_by(Trt.1, R, Trt.2)%>%
    mutate( Y = round(rnorm(n(), mu.y, sigma.star), 2) )

  
   check <-  group_by(dat, Trt.1, R, Trt.2)%>%
    summarise(check.mu = mean(Y),
              check.var = var(Y),
              N.obs = n())
  
  check.2 <- group_by(dat, Trt.1)%>%
          summarise(N.Trt.1 = n())
  

  
  check.joined <- left_join(check, check.2, by="Trt.1")

  
  check.joined <- filter(check.joined, is.na(Trt.2))%>%
                    mutate(res.rate = round(N.obs/N.Trt.1,2))

  
  # Estimated g1 and g2 from data: added 9 MAR 2018:
  
  g1.est <- check.joined$res.rate[check.joined$Trt.1 == "A"]
  g2.est <- check.joined$res.rate[check.joined$Trt.1 == "B"]
  

  g1 <- g1.est
  g2 <- g2.est

  

  simu.out <- list()

  # calling the function: 
  # to_test <- function(data, treatment.seq, placebo.seq, indi.shared.dist, indi.noninf.equiv )
  
  # for non-inferiority test:
  
  power.out.Di.NI <- to_test(dat, c("A", "C"), c("B", "C"), "Di", "NI", N, g1, g2, sigma,
                             theta, z_alpha)

  
  power.out.Sh.NI <- to_test(dat, c("B", "D"), c("B", "C"), "Sh", "NI", N, g1, g2, sigma, 
                             theta, z_alpha)


  
    
  simu.out$I.Di.NI <- power.out.Di.NI$I
  simu.out$Y.treat.avg.Di.NI <- power.out.Di.NI$Y.treat.avg
  simu.out$Y.placebo.avg.Di.NI <- power.out.Di.NI$Y.placebo.avg
  
  
  simu.out$I.Sh.NI <- power.out.Sh.NI$I
  simu.out$Y.treat.avg.Sh.NI <- power.out.Sh.NI$Y.treat.avg
  simu.out$Y.placebo.avg.Sh.NI <- power.out.Sh.NI$Y.placebo.avg

  
  # for equivalence test:
  
  power.out.Di.EQ <- to_test(dat, c("A", "C"), c("B", "C"), "Di", "EQ",  N, g1, g2, sigma,
                             theta, z_alpha)

  
  power.out.Sh.EQ <- to_test(dat, c("B", "D"), c("B", "C"), "Sh", "EQ",  N, g1, g2, sigma,
                             theta, z_alpha)


  simu.out$I.Di.EQ <- power.out.Di.EQ$I
  simu.out$Y.treat.avg.Di.EQ <- power.out.Di.EQ$Y.treat.avg
  simu.out$Y.placebo.avg.Di.EQ <- power.out.Di.EQ$Y.placebo.avg
  
  
  simu.out$I.Sh.EQ <- power.out.Sh.EQ$I
  simu.out$Y.treat.avg.Sh.EQ <- power.out.Sh.EQ$Y.treat.avg
  simu.out$Y.placebo.avg.Sh.EQ <- power.out.Sh.EQ$Y.placebo.avg
  
 
  return(simu.out)

  
}

