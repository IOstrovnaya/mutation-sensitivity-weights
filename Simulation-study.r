#________________________________________________________
#___ A Weighted Estimation Approach to Correct for 
#      Imperfect Mutation Detection Due to Intra-tumor 
#      Heterogeneity in Tumor Sequencing
#___ Audrey Mauguen, Irina Ostrovnaya
#___ October 2026
#________________________________________________________

# This code is a companion for the manuscript A Weighted Estimation Approach to Correct for Imperfect Mutation Detection 
# Due to Intra-tumor Heterogeneity in Tumor Sequencing


# The question is: how to properly estimate and compare mutations rates between two groups where the chances of finding 
#  the mutations are imperfect and not equal between the groups

source("Functions-for-simulations.r")


#___  Parameters set up

nsim <- 1000    # number of simulations
tres <- NULL    # initialization of the results table
set.seed(159)   # seed for reproducibility

# Parameters for the Beta distribution of the parameter theta (=heterogeneity, ie chances of observing a mutation if present)
par1 <- c(10, 20)    # average theta = 0.33
par2 <- c(15, 15)    # average theta = 0.50
par3 <- c(1, 0)      # average theta = 1
par.beta <- list(par1, par2, par3)


for(N in c(200, 400, 1000)) for(pm in c(0.1, 0.5)) for(pr in 0.5) for(OR in c(1, 1.5, 2)) for(th0.beta.par in par.beta) for(th1.beta.par in par.beta){
 
  #___ Initiation of the estimated parameters to be saved; 
  # or = Odds ratio; 
  # b1 and b1.se = parameter estimate (log(OR)) with standard error
  # or.nw, b1.nw and b1.se.nw: same parameters estimated with a non weighting model
  # p.wtest: p-value to test if OR != 0, from Wald test, usign the weighted model
  # p.nwtest: Wald p-value as above but from the non-weighted model
  # or.uncond, b1.uncond, b1.se.uncond and p.wtest.uncond: same parameters using the 
  #  "unconditional" mean for theta, ie correcting for the fact that theta is only measured 
  #  in patients with mutation
  or <- b1 <- b1.se <- or.nw <- b1.nw <- b1.se.nw <- p.wtest <- p.nwtest  <- rep(NA, nsim) 
  or.uncond <- b1.uncond <- b1.se.uncond <-  p.wtest.uncond <- rep(NA, nsim) 
  # Average theta parameters (=heterogeneity; ie chance of seing a mutation that is present)
  th0 <- th0.beta.par[1] / (th0.beta.par[1] + th0.beta.par[2])
  th1 <- th1.beta.par[1] / (th1.beta.par[1] + th1.beta.par[2])
  # 0 = Non-Reponders;  1 = Responders
  
    for(s in 1:nsim){
      print(paste("s=", s))
      
      #_____ Data simulation (see "Function-for-simulations.r")
       ds <- datasimOR(pm, pr = pr, OR = OR, th0=NULL, th1=NULL, N, theta = "beta", sd.theta = NULL, 
                            th0.beta.par, th1.beta.par)
        
      #_____  Analysis
       
      # requires a minimum number of mutation in each group for model estimation
      if(sum(ds$mut[ds$resp == 0]) > 3 & sum(ds$mut[ds$resp == 1]) > 3){
        
        #__ estimating the average theta
        m.theta0 <- mean(ds$heterog[ds$mut == 1 & ds$resp == 0], na.rm = TRUE)
        m.theta1 <- mean(ds$heterog[ds$mut == 1 & ds$resp == 1], na.rm = TRUE)
        
        # Weights for each subject
        ds$weighti <- NA
        ds$weighti[ds$mut == 0 & ds$resp == 0] <- m.theta0
        ds$weighti[ds$mut == 0 & ds$resp == 1] <- m.theta1
        # If subject has an observed mutation -> weight = 1
        ds$weighti[ds$mut == 1] <- 1
        
        # using unconditional mean (ie, corrected for the fact that theta is observed in 
        #  subjects with observed mutation only)
        m.theta0.uncond <- 1/(mean(1/ds$heterog[ds$mut == 1 & ds$resp == 0], na.rm = TRUE))
        m.theta1.uncond <- 1/(mean(1/ds$heterog[ds$mut == 1 & ds$resp == 1], na.rm = TRUE))
        
        # Weights for each subject, using unconditional mean
        ds$weighti.uncond <- NA
        ds$weighti.uncond[ds$mut == 0 & ds$resp == 0] <- m.theta0.uncond
        ds$weighti.uncond[ds$mut == 0 & ds$resp == 1] <- m.theta1.uncond
        # If subject has an observed mutation -> weight = 1
        ds$weighti.uncond[ds$mut == 1] <- 1
        
        if(!is.na(sum(ds$weighti.uncond))){ # check that the weight is available
        
          # Estimation of the effect using weighted logistic regression
          lrw <- glm(resp ~ mut, data = ds, weights = weighti, family = binomial(link = "logit"))
          # Estimation of the effect using the weighted logistic regression, using the unconditional mean of theta
          lrw.uncond <- glm(resp ~ mut, data = ds, weights = weighti.uncond, family = binomial(link = "logit"))

          # Save the estimates in a vector
          or[s] <- round(exp(lrw$coefficients[2]), 3)
          b1[s] <- round(lrw$coefficients[2], 3)
          b1.se[s] <- summary(lrw)$coeff[2, 2]

          or.uncond[s] <- round(exp(lrw.uncond$coefficients[2]), 3)
          b1.uncond[s] <- round(lrw.uncond$coefficients[2], 3)
          b1.se.uncond[s] <- summary(lrw.uncond)$coeff[2, 2]
          
          p.wtest[s]  <- summary(lrw)$coeff[2, 4]
          p.wtest.uncond[s]  <- summary(lrw.uncond)$coeff[2, 4]
        
          # Estimation of the effect using a non-weighted logistic regression
          lr <- glm(resp ~ mut, data = ds, family = binomial)
          
          # Save the estimates in a vector
          or.nw[s] <- round(exp(lr$coefficients[2]), 3)
          b1.nw[s] <- round(lr$coefficients[2], 3)
          b1.se.nw[s] <- summary(lr)$coeff[2, 2]
          
          p.nwtest[s]  <- summary(lr)$coeff[2, 4]
        
         }
      } }
    # Save the results in an overall table
    tres <- rbind(tres,
                  cbind(# parameters for each scenario
                        pm = pm, pr = pr, OR = OR, beta1 = round(log(OR), 3),
                        th0 = th0, th1 = th1,
                        N = N, 
                        # average estimates for the scenario, with power
                        or = median(or, na.rm = T), b1 = mean(b1, na.rm = T), b1.se = mean(b1.se, na.rm = T), 
                        power = sum(p.wtest <= 0.05, na.rm=T) / sum(!is.na(p.wtest)), 
                        or.uncond = median(or.uncond, na.rm = T), b1.uncond = mean(b1.uncond, na.rm = T), b1.se.uncond = mean(b1.se.uncond, na.rm = T), 
                        power.uncond = sum(p.wtest.uncond <= 0.05, na.rm=T) / sum(!is.na(p.wtest.uncond)), 
                        or.nw = median(or.nw, na.rm = T), b1.nw = mean(b1.nw, na.rm = T), b1.se.nw = mean(b1.se.nw, na.rm = T), 
                        power.nw = sum(p.nwtest <= 0.05, na.rm=T) / sum(!is.na(p.nwtest)),
                        # Number of times the model did not converge (number of NAs)
                        n.na = sum(is.na(or)), n.na.uncond = sum(is.na(or.uncond)), n.na.nw = sum(is.na(or.nw)) 
                        ))
}


save(tres, file = "tres.Rdata")

