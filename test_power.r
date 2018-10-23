  


  # function to test both non-inferiority and superiority:

  to_test <- function(data, treatment.seq, placebo.seq, indi.shared.dist, indi.noninf.equiv, N, g1, g2, sigma, 
                      theta, z_alpha){

    
    # indi.shared.dist: "Sh" for shared and "Di" for distint path
    # indi.noninf.equiv: "NI" for non-inferiority and "EQ" for equivalence test
  
  
  ## responder: higher the better
  # responder and non-responder with response rate g1 and g2 in two gropus:
    
  Y.treat.non_res = filter(data, (Trt.1 == treatment.seq[1] & Trt.2 == treatment.seq[2]))%>%
    ungroup()%>%select(Y)%>%unlist(use.names = FALSE)
  
  Y.treat.res = filter(data, (Trt.1 == treatment.seq[1] & is.na(Trt.2) == T))%>%
    ungroup()%>%select(Y)%>%unlist(use.names = FALSE)
  
  Y.treat <- c(Y.treat.res, Y.treat.non_res)
  
  Y.placebo.non_res = filter(data, (Trt.1 == placebo.seq[1] & Trt.2 == placebo.seq[2]))%>%
    ungroup()%>%select(Y)%>%unlist(use.names = FALSE)
  Y.placebo.res =  filter(data, (Trt.1 == placebo.seq[1] & is.na(Trt.2) == T))%>%
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
  
  mu.max <- max(mean(Y.treat.res), mean(Y.treat.non_res), mean(Y.placebo.res), mean(Y.placebo.non_res))
  mu.min <- min(mean(Y.treat.res), mean(Y.treat.non_res), mean(Y.placebo.res), mean(Y.placebo.non_res))
  
  
  # estimate sigma: 
  
  n.treat <- length(Y.treat)
  n.placebo <- length(Y.placebo)
  sigma.1 <- sqrt(( (n.treat-1)*var(Y.treat) + (n.placebo -1)*var(Y.placebo))/(n.treat + n.placebo -2))
  
  if(is.na(sigma.1)){ print(c("sigma.1", sigma.1))}
  
  if(is.na(sigma.1)){ sigma <- sigma}
  if(!is.na(sigma.1)){ sigma <- sigma.1}
  
  est.mu.AA <- mean(Y.treat.res)
  est.mu.AC <- mean(Y.treat.non_res)
  
  est.mu.BB <- mean(Y.placebo.res)
  est.mu.BC <- mean(Y.placebo.non_res)
  
  # for shared path only:
  
  est.mu.BD <- mean(Y.treat.non_res)
  
  

  if(indi.shared.dist == "Di"){  
    

    V = ( 2*(4 - g1 -g2)*(1/N)*(sigma^2) + (1/N)*(  g1*(2-g1)*(est.mu.AA^2) + (3-2*g1-g1^2)*(est.mu.AC^2) -2*g1*(1-g1)*est.mu.AA*est.mu.AC
                                                  + g2*(2-g2)*(est.mu.BB^2) + (3-2*g2-g2^2)*(est.mu.BC^2) -2*g2*(1-g2)*est.mu.BB*est.mu.BC ) )
    
    if(is.na(V)){ print(c("V is NA", V))}
    if(is.infinite(V)){ print(c("V is infinite", V))}
    if(is.nan(V)){ print(c("V is nan", V))}

  }
  
  
  if(indi.shared.dist == "Sh"){  
    
    if(treatment.seq[1] == "B" & placebo.seq[1] == "B"){
      

      V = ( 2*(4 - g2 -g2)*(1/N)*(sigma^2) + (1/N)*(  g2*(2-g2)*(est.mu.BB^2) + (3-2*g2-g2^2)*(est.mu.BC^2) -2*g2*(1-g2)*est.mu.BB*est.mu.BC
                                                    + g2*(2-g2)*(est.mu.BB^2) + (3-2*g2-g2^2)*(est.mu.BD^2) -2*g2*(1-g2)*est.mu.BB*est.mu.BD ) 
            - 2*(1/N)*(2*g2*(sigma^2 + est.mu.BB^2) - (g2*est.mu.BB + (1-g2)*est.mu.BC)*(g2*est.mu.BB + (1-g2)*est.mu.BD)) )
      
      
      if(is.na(V)){ print(c("V is NA", V))}
      if(is.infinite(V)){ print(c("V is infinite", V))}
      if(is.nan(V)){ print(c("V is nan", V))}
      if(V <= 0){ print(c("V is 0 or negetive", V))}

      }
  }
  
  
  # for non-inferiority test:
  
  if(indi.noninf.equiv == "NI"){
    Z.non_inf = (Y.treat.avg - Y.placebo.avg + theta)/(sqrt(V) )  
        
    I.non_inf = 1*(Z.non_inf >   ( z_alpha  ))
    I = I.non_inf
    
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
    

    
    if(is.na(I) | is.na(I.non_inf) | is.na(I.equiv) | is.na(V) | is.na(Y.treat.avg) | is.na(Y.placebo.avg) ){
      print("NA in Equivalence")
      print(c(I, I.non_inf, I.equiv, V, Y.treat.avg, Y.placebo.avg)) }
    
  }
  
  simu.out <- list()
  simu.out$I <- I
  simu.out$Y.treat.avg <- Y.treat.avg
  simu.out$Y.placebo.avg <- Y.placebo.avg
  
  return(simu.out)
  
  } # end of to_test function
  

  