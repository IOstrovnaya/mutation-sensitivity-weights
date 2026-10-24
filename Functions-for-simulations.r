#________________________________________________________
#___ A Weighted Estimation Approach to Correct for 
#      Imperfect Mutation Detection Due to Intra-tumor 
#      Heterogeneity in Tumor Sequencing
#___ Audrey Mauguen, Irina Ostrovnaya
#___ October 2026
#________________________________________________________


#___ Functions to simulate the data

# pm:   probability of mutation overall
# pr:   probability of response/being a case (ie, probability of being in group 1)
# OR:   Odds Ratio, effect of the mutation on the outcome (case/control or responder/non-responder)
# th0:  heterogeneity in group 0 (controls, non responders), ie probability to observe the mutation if present
# th1:  heterogeneity in group 1 (cases, responders), ie probability to observe the mutation if present
# 
# N:   number of subjects

# theta: distribution of theta within the two groups
#    - "same": all subjects have the same fixed theta, taking the value th0 or th1
#    - "normal": theta follows a Normal distribution, with average th0 or th1
#      -- sd.theta: standard-deviation for the distribution of theta, identical in both groups
#    - "beta": theta follows a Beta distribution, with shape parameters given in arguments 
#      -- th0.beta.par: two shape parameters for the Beta distribution in group 0
#      -- th1.beta.par: two shape parameters for the Beta distribution in group 1

datasimOR <- function(pm, pr, OR = 1, th0, th1, N = 200, theta = "same", sd.theta = 0.25, 
                      th0.beta.par = c(NULL, NULL), th1.beta.par = c(NULL, NULL)){
  # Subject id
  pt <- 1:N
  
  # mutation -- true status
  m <- rbinom(N, 1, pm) #patient level mutation status
  N.m1 <- sum(m)        # number of mutate patients
  N.m0 <- N - N.m1      # number of non-mutated patients
  
  # responses - using logistic model
  b0 <- -log((1 - pr) / pr)  # average response - intercept
  b1 <- log(OR)              # average response - mutation impact
  pri <- 1 / (1 + exp(-(b0 + b1*m)) ) # individual probability of response based on the mutation status
  resp <- rbinom(N, 1, pri)  # individual response status based on the probability of reponding
  
  N.r1 <- sum(resp)   # number of responders
  N.r0 <- N - N.r1    # number of non-responders
  
  # Simulation of theta
  heterog <- NA
  
  if(theta == "same"){
    heterog[resp == 1] <- rep(th1, N.r1)
    heterog[resp == 0] <- rep(th0, N.r0)
  }
  
  if(theta == "normal"){
    heterog[resp == 1] <- rnorm(n = N.r1, mean = th1, sd = sd.theta)
    heterog[resp == 0] <- rnorm(n = N.r0, mean = th0, sd = sd.theta)
    #_ Truncate values not between 0 and 1
    heterog[heterog < 0] <- 0.001; heterog[heterog > 1] <- 1
  }
  
  if(theta == "beta"){
    heterog[resp == 1] <- rbeta(n = N.r1, shape1 = th1.beta.par[1], shape2 = th1.beta.par[2])
    heterog[resp == 0] <- rbeta(n = N.r0, shape1 = th0.beta.par[1], shape2 = th0.beta.par[2])
    
    th0 <- th0.beta.par[1] / (th0.beta.par[1] + th0.beta.par[2])
    th1 <- th1.beta.par[1] / (th1.beta.par[1] + th1.beta.par[2])
  }
  
  # Mutation -- observed status
  #  mutation present X mutation observed
  mut <- NA
  mut[resp == 1] <- m[resp == 1] * (runif(N.r1) <= th1)
  mut[resp == 0] <- m[resp == 0] * (runif(N.r0) <= th0)
  
  d <- as.data.frame(cbind(pt = pt, mut = mut,  resp = resp, heterog = heterog, truem = m))
  
  return(d)
}






