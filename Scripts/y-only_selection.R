##### Modelling the y-only process with selection#####



#Choose a parent with y1, y2 many copies of an allele according to the current distribution of y1,y2
#With probability r recombination may produce a new number of copies, otherwise either y1 or y2
#The recombination is modelled as the sum of two independent uniform(y1,y2) distributed random variables (break points on the chromosomes)
#Iterating this process yields to a stationary distribution of y



#####Functions#####

#Selection
approx_sel <- function(y,sx,sy,eps){
  #Quadratic Taylor approximation
  s_tilde <- eps*exp(-eps)*sqrt(sx*sy)
  y_opt <- ((log(sx)-log(sy))/(2*eps)+1)
  out <- max(0,1-s_tilde*((y-y_opt)^2))
  return(out)
}


#Set of possible copy numbers
y_set <- function(sx,sy,eps){
  s_tilde <- eps*exp(-eps)*sqrt(sx*sy)
  y_opt <- ((log(sx)-log(sy))/(2*eps)+1)
  #Calculate the set of copy numbers with positive fitness
  min_y <- max(1,as.integer(y_opt-1/sqrt(s_tilde))+1)
  max_y <- as.integer(y_opt+1/sqrt(s_tilde))
  out <- min_y:max_y
  return(out)
}


#Trapez function, as the probability to produce a gamete with k gene copies from parental chromosomes with y1, y2 genes
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
  n <- length(k)
  x_tmp <- rep(0,n)
  
  #f_post(k) = SUM_{y1,y2} f_prior(y1)*f_prior(y2)*P[k | y1,y2]
  #where P is the probability to create a gamete with k many copies out of the parental chromosomes with y1,y2 many copies, defined in function conditional_prob

  Z <- sum(f_prior*approx_sel(k,sx,sy,eps))
  f_sel <- f_prior*approx_sel(k,sx,sy,eps)/Z #normalisation
  
  for (y1 in 1:n) {
    for(y2 in 1:n) {
      x_tmp <- x_tmp + f_sel[y1]*f_sel[y2]*conditional_prob(k,r,y1,y2) 
    }
    
  }
  return(x_tmp)
}








###### Initial Parameters #####

sx <- 0.02
sy <- 0.005
eps <- 0.05
n <- length(y_set(sx,sy,eps))
r <- 0.01

#blank plot
plot(y_set(sx,sy,eps),1:n/n,type="n",ylim = c(0,0.25),xlab="Copy number",ylab = "Frequency")


#run 3 different scenarios
for (i in 1:3) {
  r <- c(0.01,0.1,0.5)[i]
  c <- c("black","blue","red")[i]

  #initial distribution
  y_opt <- 1/2*((log(sx)-log(sy))/(2*eps)+1)
  current_distr <- dnorm(y_set(sx,sy,eps),mean=y_opt,sd=y_opt/3)

  f_tmp <- f_post(y_set(sx,sy,eps),0.1,current_distr)
  difference <- 1
  
  #iterate until convergence
  while (difference>0.001) {
  f_tmp <- f_post(y_set(sx,sy,eps),r,current_distr) #update by function f_post (posterior)
  difference <- sum(abs(1/sum(f_tmp)*f_tmp-current_distr))
  current_distr <- 1/sum(f_tmp)*f_tmp #normalise
}



  #plot the stationary distribution
  points(y_set(sx,sy,eps),current_distr,type="l",col=c)

  #print the mean value of the distribution for a given recombination rate
  print(c(r,sum(y_set(sx,sy,eps)*current_distr)))
}

#mark the calculated optimal value
abline(v=y_opt,col="red")

y_opt