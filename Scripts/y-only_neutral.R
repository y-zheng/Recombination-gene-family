###### Modelling the y-only process without selection #####



#Choose a parent with y1, y2 many copies of an allele according to the current distribution of y1,y2
#With probability r recombination may produce a new number of copies, otherwise either y1 or y2
#The recombination is modelled as the sum of two independent uniform(y1,y2) distributed random variables (break points on the chromosomes)
#Iterating this process yields to a stationary distribution of y, depending on the Expectation of the initial distribution
#Increasing r speeds up the time to stationary distribution, but does not affect it



#####Functions#####

#Trapez function
trapez <- function(k,y1,y2){
  x_tmp <- rep(0,length(k))
  #Choose one break point on each chromosome, i.e. any of the alleles
  #By uneven recombination, the outcome is randomly distributed on 1,...,y1+y2-1 according to the density of the sum of two uniform distributed random variables, which is the Trapez-distribution (for y1=y2 it is the Triangle distribution)
  x_tmp[k<=y1 & k<=y2] <- k[k<=y1 & k<=y2]
  x_tmp[k<=y1 & k>y2] <- y2
  x_tmp[k>y1 & k<=y2] <- y1
  x_tmp[k>y1 & k>y2] <- y1+y2-k[k>y1 & k>y2]
  x_tmp[k<1] <- 0
  x_tmp[k>(y1+y2-1)] <- 0
  
  return(x_tmp/(y1*y2))
}

#Posterior function

conditional_prob <- function(k,r,y1,y2){
  x_tmp <- rep(0,length(k))
  
  #Given parent with y1,y2 many copies, the gamete is with probability 1/2 y1,y2 if no recombination happens
  x_tmp[k==y1] <- 1/2*(1-r)
  x_tmp[k==y2] <- 1/2*(1-r)
  x_tmp[k==y1 & k==y2] <- (1-r)
  #If recombination happens, the probability of k as an outcome for given y1,y2 is the sum of two uniformly distributed random variables, which is trapez-distributed (trapez function defined down below)
  x_tmp <- x_tmp + r*trapez(k,y1,y2)
  
  return(x_tmp)
}

f_post <- function(k,r,f_prior) {#return f_post(k), for given recombination rate r and prior distribution f_prior
  x_tmp <- rep(0,length(k))
  #f_post(k) = SUM_{y1,y2} f_prior(y1)*f_prior(y2)*P[k | y1,y2]
  #where P is the probability to create a gamete with k many copies out of the parental chromosomes with y1,y2 many copies, defined in function conditional_prob
  for (y1 in 1:length(f_prior)) {
    for(y2 in 1:length(f_prior)) {
      x_tmp <- x_tmp + f_prior[y1]*f_prior[y2]*conditional_prob(k,r,y1,y2) 
    }
    
  }
  return(x_tmp)
}








###### Initial Parameters #####

n <- 100 #max length (should be infinity, but we ignore the possibility to create more than n copies by recombination)
r <- 0.01 #recombination rate


#initial distribution:
current_distr <- rep(0,n) 

current_distr[2] <- 1/2
current_distr[18] <- 1/2

current_distr <- dnorm(1:n,mean=10,sd=3)/sum(dnorm(1:n,mean=10,sd=3))


#stationary distribution:
E0 <- sum(1:n * current_distr)
stat_distr <- (1:n)*exp(-2/E0*1:n)



##### Iterating #####

plot(1:n,current_distr,type="l",main="Copy number Distribution",xlab = "number of copies",ylab="Frequency",xlim=c(0,40),ylim=c(0,0.15))


for(i in 1:100) { #iterate 100x
  f_tmp <- f_post(1:n,r,current_distr) #update by function f_post (posterior)
  current_distr <- 1/sum(f_tmp)*f_tmp #normalise (since we ignore the values >n)
}

points(1:n,current_distr,type="l",col="blue",lty=3)



for(i in 1:150) { #iterate 150x
  f_tmp <- f_post(1:n,r,current_distr) #update by function f_post (posterior)
  current_distr <- 1/sum(f_tmp)*f_tmp #normalise (since we ignore the values >n)
}

points(1:n,current_distr,type="l",col="green",lty=2)


for(i in 1:750) { #iterate 750x
  f_tmp <- f_post(1:n,r,current_distr) #update by function f_post (posterior)
  current_distr <- 1/sum(f_tmp)*f_tmp #normalise (since we ignore the values >n)
}


points(1:n,current_distr,type="l",col="red")

points(1:n,stat_distr/(sum(stat_distr)))
