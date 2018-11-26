  


  # function to test both non-inferiority and superiority:



  to_test <- function(data, treatment.seq, placebo.seq, indi.shared.dist, indi.noninf.equiv, N, g.treat, g.placebo,
                      theta, z_alpha){

    
    # indi.shared.dist: "Sh" for shared and "Di" for distint path
    # indi.noninf.equiv: "NI" for non-inferiority and "EQ" for equivalence test
  
  
  ## responder: higher the better
  # responder and non-responder with response rate g.treat and g.placebo in two gropus:
    
  Y.treat.non_res = filter(data, (Trt.1 == treatment.seq[1] & Trt.2 == treatment.seq[2]))%>%
    ungroup()%>%select(Y)%>%unlist(use.names = FALSE)
  
  Y.treat.res = filter(data, (Trt.1 == treatment.seq[1] & R == 1))%>%
    ungroup()%>%select(Y)%>%unlist(use.names = FALSE)
  
  Y.treat <- c(Y.treat.res, Y.treat.non_res)
  
  Y.placebo.non_res = filter(data, (Trt.1 == placebo.seq[1] & Trt.2 == placebo.seq[2]))%>%
    ungroup()%>%select(Y)%>%unlist(use.names = FALSE)
  Y.placebo.res =  filter(data, (Trt.1 == placebo.seq[1] & R == 1))%>%
    ungroup()%>%select(Y)%>%unlist(use.names = FALSE)
  Y.placebo <- c(Y.placebo.res, Y.placebo.non_res) 
  
  

  
  if(is.nan(mean(Y.treat.res)) == T ){Y.treat.res <- 0}
  if(is.nan(mean(Y.placebo.res)) == T ){Y.placebo.res <- 0}
  
  if(is.nan(mean(Y.treat.non_res)) == T ){Y.treat.non_res <- 0}
  if(is.nan(mean(Y.placebo.non_res)) == T ){Y.placebo.non_res <- 0}
  
  if(is.infinite(mean(Y.treat.res)) == T ){Y.treat.res <- 0}
  if(is.infinite(mean(Y.placebo.res)) == T ){Y.placebo.res <- 0}
  
  if(is.infinite(mean(Y.treat.non_res)) == T ){Y.treat.non_res <- 0}
  if(is.infinite(mean(Y.placebo.non_res)) == T ){Y.placebo.non_res <- 0}
  
  ## group wise weighted mean outcome:
  
  Y.treat.avg = (1/N)*(4*sum(Y.treat.non_res) + 2*sum(Y.treat.res) )
  Y.placebo.avg = (1/N)*(4*sum(Y.placebo.non_res) + 2*sum(Y.placebo.res) )
  
  if( is.infinite(Y.placebo.avg)){
    print("NA in Y.treat.avg, Y.placebo.avg")
    print(c(Y.treat.avg, Y.placebo.avg))  }
  

  # estimate sigma: 
  
  n.treat <- length(Y.treat)
  n.placebo <- length(Y.placebo)
  sigma.1 <- sqrt(( (n.treat-1)*var(Y.treat) + (n.placebo -1)*var(Y.placebo))/(n.treat + n.placebo -2))
  
  if(is.na(sigma.1)){ print(c("sigma.1", sigma.1))}
  
  if(is.na(sigma.1)){ 
    # sigma <- sigma
    print("Sigma is NA!")
    }
  if(!is.na(sigma.1)){ sigma <- sigma.1}
  
  est.mu.treat.res <- mean(Y.treat.res)
  est.mu.treat.non_res <- mean(Y.treat.non_res)
  
  est.mu.placebo.res <- mean(Y.placebo.res)
  est.mu.placebo.non_res <- mean(Y.placebo.non_res)
  


  if(indi.shared.dist == "Di"){  
    

    V = ( 2*(4 - g.treat -g.placebo)*(1/N)*(sigma^2) + (1/N)*(  g.treat*(2-g.treat)*(est.mu.treat.res^2) + (3-2*g.treat-g.treat^2)*(est.mu.treat.non_res^2) -2*g.treat*(1-g.treat)*est.mu.treat.res*est.mu.treat.non_res
                                                  + g.placebo*(2-g.placebo)*(est.mu.placebo.res^2) + (3-2*g.placebo-g.placebo^2)*(est.mu.placebo.non_res^2) -2*g.placebo*(1-g.placebo)*est.mu.placebo.res*est.mu.placebo.non_res ) )
    
    if(is.na(V)){ print(c("V is NA", V))}
    if(is.infinite(V)){ print(c("V is infinite", V))}
    if(is.nan(V)){ print(c("V is nan", V))}

  }
  
  
  if(indi.shared.dist == "Sh"){  
    


      V = ( 2*(4 - g.placebo -g.placebo)*(1/N)*(sigma^2) + (1/N)*(  g.placebo*(2-g.placebo)*(est.mu.placebo.res^2) + (3-2*g.placebo-g.placebo^2)*(est.mu.placebo.non_res^2) -2*g.placebo*(1-g.placebo)*est.mu.placebo.res*est.mu.placebo.non_res
                                                    + g.placebo*(2-g.placebo)*(est.mu.placebo.res^2) + (3-2*g.placebo-g.placebo^2)*(est.mu.treat.non_res^2) -2*g.placebo*(1-g.placebo)*est.mu.placebo.res*est.mu.treat.non_res ) 
            - 2*(1/N)*(2*g.placebo*(sigma^2 + est.mu.placebo.res^2) - (g.placebo*est.mu.placebo.res + (1-g.placebo)*est.mu.placebo.non_res)*(g.placebo*est.mu.placebo.res + (1-g.placebo)*est.mu.treat.non_res)) )
      
      
      if(is.na(V)){ print(c("V is NA", V))}
      if(is.infinite(V)){ print(c("V is infinite", V))}
      if(is.nan(V)){ print(c("V is nan", V))}
      if(V <= 0){ print(c("V is 0 or negetive", V))}

      # }
  }
  
  
  # for non-inferiority test:
  
  if(indi.noninf.equiv == "NI"){
    Z.non_inf = (Y.treat.avg - Y.placebo.avg + theta)/(sqrt(V) )  
        
    I.non_inf = 1*(Z.non_inf >   ( z_alpha  ))
    I = I.non_inf
    
    if(I == 1){ Decision = "Reject the null hypothesis" 
    } else { Decision = "Unable to reject the null hypothesis" }
    
    
    
    p.val.Non_Inferiority <- round(pnorm(-abs(Z.non_inf)), 5)
    upp.bound.Bayes.fac.NI <- 1/(-exp(1)*p.val.Non_Inferiority*log(p.val.Non_Inferiority))
    
    p.val.Superiority <- "-"
    upp.bound.Bayes.fac.Sup <- "-"
    
    if(is.na(I) | is.na(I.non_inf) | is.na(V) | is.na(Y.treat.avg) | is.na(Y.placebo.avg) ){
      print("NA in Non-inferiority")
    print(c(I, I.non_inf, V, Y.treat.avg, Y.placebo.avg)) }
    }
  
  # for equivalence test:
  
  if(indi.noninf.equiv == "EQ"){
    Z.non_inf = (Y.treat.avg - Y.placebo.avg + theta)/(sqrt(V) )  
    Z.equiv   = (Y.treat.avg - Y.placebo.avg - theta)/(sqrt(V) )  
    
    I.non_inf = 1*(Z.non_inf >   ( z_alpha  ))
    I.equiv   = 1*(Z.equiv   <   (-z_alpha  ))
    
    I = I.non_inf*I.equiv
    
    if(I == 1){ Decision = "Reject the null hypothesis" 
    } else { Decision = "Unable to reject the null hypothesis" }
    
    
    

    
    p.val.Non_Inferiority <- round(pnorm(-abs(Z.non_inf)), 5)
    upp.bound.Bayes.fac.NI <- 1/(-exp(1)*p.val.Non_Inferiority*log(p.val.Non_Inferiority))
    
    p.val.Superiority <- round(pnorm(-abs(Z.equiv)), 5)
    upp.bound.Bayes.fac.Sup <- 1/(-exp(1)*p.val.Superiority*log(p.val.Superiority))
    

    
    if(is.na(I) | is.na(I.non_inf) | is.na(I.equiv) | is.na(V) | is.na(Y.treat.avg) | is.na(Y.placebo.avg) ){
      print("NA in Equivalence")
      print(c(I, I.non_inf, I.equiv, V, Y.treat.avg, Y.placebo.avg)) }
    
  }
  
  simu.out <- list()
  simu.out$NI_EQ <- indi.noninf.equiv
  simu.out$distint.shared.path <- indi.shared.dist
  # simu.out$I <- I
  simu.out$Y.treat.avg <- round(Y.treat.avg, 2)
  simu.out$Y.placebo.avg <- round(Y.placebo.avg, 2)
  # simu.out$p.val <- p.val
  simu.out$p.val.Non_Inferiority <- p.val.Non_Inferiority
  simu.out$upp.bound.Bayes.fac.NI <- ifelse(is.na(upp.bound.Bayes.fac.NI), Inf, upp.bound.Bayes.fac.NI)
  
  simu.out$p.val.Superiority <- p.val.Superiority
  simu.out$upp.bound.Bayes.fac.Sup <- ifelse(is.na(upp.bound.Bayes.fac.Sup), Inf, upp.bound.Bayes.fac.Sup)
  
  simu.out$Decision <- Decision
  
  # print(data.frame(unlist(simu.out)))
  
  
  return(simu.out)
  
  } # end of to_test function
  

  