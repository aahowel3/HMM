#the binomial equation itself here is fine
#prob of being in either/or component is a binary 0 or 1 thats still good - until we try it with 3 components 
#but what needs to change is the likelihood dist our samples are coming from
#not 2 binomial dists governed by a single param switching based on a binomial eq
#but 2 normal dists goverened by 2 params swithching based on a binomial eq 
expectation <- function(sample,p,a1,a2,b1,b2)
{
  p_expectation <- (p*dnorm(sample,m=a1,s=b1)) / ( (p*dnorm(sample,m=a1,s=b1)) + ((1-p)*dnorm(sample,m=a2,s=b2)) )
  return(p_expectation)
}


maximization <- function(sample,epart){
  
  # estimate p
  
  p_temp <- mean(epart)
  
  # estimate a and b
#shit this needs to be different - estimating 2as and 2bs  
  #mu is also a different calculation from sigma 
#http://sia.webpopix.org/mixtureModels.html#the-em-algorithm   
#these calculations will probably change when using a negbionmial    
  
  a1_temp <- sum(sample*epart) / sum(epart)
  b1_temp = sqrt(sum(epart*(sample^2))/sum(epart)-(a1_temp)^2)
  
  a2_temp <- sum(sample*(1-epart)) / sum(1-epart)
  epart2=1-epart
  b2_temp = sqrt(sum(epart2*(sample^2))/sum(epart2)-(a2_temp)^2)
  
  
  list(p_temp,a1_temp,a2_temp,b1_temp,b2_temp)   
}




EM <- function(sample,p_init,a1_init,a2_init,b1_init,b2_init,maxit=1000,tol=1e-6)
{
  # Estimation of parameter(Initial)
  flag <- 0
  p_cur <- p_init; a1_cur <- a1_init; a2_cur <- a2_init; b1_cur <- b1_init; b2_cur <- b2_init
  
  # Iterate between expectation and maximization parts
  
  for(i in 1:maxit){
    cur <- c(p_cur,a1_cur,a2_cur,b1_cur,b2_cur)
    new <- maximization(sample,expectation(sample, p_cur, a1_cur,a2_cur,b1_cur,b2_cur))
  #max output order is list(p_temp,a1_temp,b1_temp,a2_temp,b2_temp)   
  
    p_new <- new[[1]]; a1_new <- new[[2]]; a2_new <- new[[3]]; b1_new <- new[[4]]; b2_new <- new[[5]]
    new_step <- c(p_new,a1_new,a2_new,b1_new,b2_new)
    
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    if( all(abs(cur - new_step) < tol) ){ flag <- 1; break}
    
    
    # Otherwise continue iteration
    p_cur <- p_new; a1_cur <- a1_new; b1_cur <- b1_new; a2_cur <- a2_new; b2_cur <- b2_new
  }
  if(!flag) warning("Didn't converge\n")
  
  list(p_cur, a1_cur, a2_cur, b1_cur, b2_cur)
}




#try as 2 normals first - see if it works 
#transition rate (p) is 0.9 (900/1000 total), a1=100, a2=5, b1=10, b2=1
MAC=rnorm(900,m=100,s=10)
IES=rnorm(100,m=5,s=1)
X = c(MAC,IES)

# Set parameter estimates
#purposefully a bit off to see if it can get closer to actuals 
p_init = 0.70; a1_init = 90; a2_init = 10; b1_init = 12; b2_init = 5


EM(X,p_init,a1_init,a2_init,b1_init,b2_init)


