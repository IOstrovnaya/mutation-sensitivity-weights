

# Data simulation based on OR
datasimOR <- function(pm, pr, OR = 1, th0, th1, N = 200, theta = "same", sd.theta = 0.25, 
                      th0.beta.par = c(NULL, NULL), th1.beta.par = c(NULL, NULL)){
  #  pm = 0.8; pr = 0.3; OR = 2; N = 200; theta = "beta"
  #  th0.beta.par = c(1, 0); th1.beta.par = c(1, 0)
  
  pt <- 1:N
  
  # mutations
  m <- rbinom(N, 1, pm) #patient level mutation status
  N.m1 <- sum(m)      # number of mutate patients
  N.m0 <- N - N.m1    # number of non-mutated patients
  
  # responses - using model
  b0 <- -log((1 - pr) / pr)
  b1 <- log(OR)
  # print(paste("OR = ", OR, " and b1 = ", b1))
  pri <- 1 / (1 + exp(-(b0 + b1*m)) )
  resp <- rbinom(N, 1, pri)
  
  N.r1 <- sum(resp)   # number of responders
  N.r0 <- N - N.r1    # number of non-responders
  
  # theta = heterogeneity
  heterog <- NA
  
  if(theta == "same"){
    heterog[resp == 1] <- rep(th1, N.r1)
    heterog[resp == 0] <- rep(th0, N.r0)
  }
  
  if(theta == "normal"){
    heterog[resp == 1] <- rnorm(n = N.r1, mean = th1, sd = sd.theta)
    heterog[resp == 0] <- rnorm(n = N.r0, mean = th0, sd = sd.theta)
    #_ Force between 0 and 1
    heterog[heterog < 0] <- 0.001; heterog[heterog > 1] <- 1
  }
  
  if(theta == "beta"){
    heterog[resp == 1] <- rbeta(n = N.r1, shape1 = th1.beta.par[1], shape2 = th1.beta.par[2])
    heterog[resp == 0] <- rbeta(n = N.r0, shape1 = th0.beta.par[1], shape2 = th0.beta.par[2])
    
    th0 <- th0.beta.par[1] / (th0.beta.par[1] + th0.beta.par[2])
    th1 <- th1.beta.par[1] / (th1.beta.par[1] + th1.beta.par[2])
  }
  mut <- NA
  mut[resp == 1] <- m[resp == 1] * (runif(N.r1) <= th1)
  mut[resp == 0] <- m[resp == 0] * (runif(N.r0) <= th0)
  
  d <- as.data.frame(cbind(pt = pt, mut = mut,  resp = resp, heterog = heterog, truem = m))
  
  return(d)
  
}

