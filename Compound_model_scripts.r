### This file contains the scripts and functions for the simulation of the "compound model" with variable copy number and allele identities.
### Note that this is a collection of functions/scripts. Do not run the whole file.
### To the end there is a few examples of ready-to-use scripts but they have to be used with the functions.

### Initialization of population size and selection, required for everything below

{
arg = commandArgs(TRUE) #Read population size and selection strength (as a multiplier to the default number of s_x = 0.01, s_y = 0.025) from command line
ne = as.numeric(arg[1])
sel_str = as.numeric(arg[2])

s1 = 0.01*sel_str
s2 = 0.0025*sel_str
k1 = 0.95
k2 = 1.05
fitlist = matrix(0,100,100) #Fitness is recorded as a matrix; copy numbers over 100 are not used because it is effectively impossible with our numbers.
fitlist[1,2] = (1+s1)
fitlist[2,2] = (1+s1)^(1+k1)
for (x in 1:100) {
 for (y in 3:100) {
  if (x <= y) {
   t1 = (1+s1)^(sum(k1^(0:(x-1))))
   t2 = (1-s2)^(sum(k2^(0:(y-3))))
   fitlist[x,y] = t1*t2
  }
 }
}

}

### Initialization of a population (blank state)

pop = list(copies = matrix(1,2*ne,5), cpn = rep(5,2*ne), fit = rep(1,ne), x = rep(1,ne), y = rep(10,ne))

### Initialization of a population (read from output)

mhc.read = function (fname,fitlist) {
 pop.o = as.matrix(read.table(fname,header=FALSE))
 cpn.o = apply(pop.o>0,1,sum)
 ne.n = nrow(pop.o)/2
 
 x = rep(0,ne.n)
 y = rep(0,ne.n)
 fit = rep(0,ne.n)
 for (i in 1:ne.n) {
  temp = c(pop.o[i,],pop.o[ne.n+i,])
  y[i] = sum(temp>0)
  x[i] = sum(unique(temp)>0)
  if (y[i] > ncol(fitlist)) {fit[i] = 0}
  else {fit[i] = fitlist[x[i],y[i]]}
 }
  
 return(list(copies = pop.o, cpn = cpn.o, fit = fit, x = x, y = y))
}

### Function to evolve a population by ONE generation

mhc.evolve = function (pop, fitlist, rec, mu, ne.n) { #Function arguments: pop = input population, fitlist = fitness matrix, rec = recombination rate, mu = mutation rate, ne.n = population size AFTER the generation
 ne = length(pop$cpn)/2
 par1 = sample(1:ne,ne.n,replace=TRUE,prob=pop$fit) #Choose parent individuals (diploid)
 par2 = sample(1:ne,ne.n,replace=TRUE,prob=pop$fit) 
 
 rec.help = c((1-rec)/2,(1-rec)/2,rec/2,rec/2)
 rec.type = sample(1:4,2*ne.n,replace=TRUE,prob=rec.help) # 1-> Gamete from paternal copy, 2-> Gamete from maternal copy, 3-> Paternal first half, maternal second half, 4-> vv
 
 pop.p = pop$copies[c(par1,par2),] #Paternal copy in father, then paternal copy in mother (for meiotic recombination)
 pop.m = pop$copies[c(par1+ne,par2+ne),]
 cpn.p = pop$cpn[c(par1,par2)]
 cpn.m = pop$cpn[c(par1+ne,par2+ne)]
 
 pop.o = matrix(0,2*ne.n,ncol(pop$copies))
 cpn.o = rep(0,2*ne.n)
 pop.o[rec.type==1,] = pop.p[rec.type==1,]
 pop.o[rec.type==2,] = pop.m[rec.type==2,]
 cpn.o[rec.type==1] = cpn.p[rec.type==1]
 cpn.o[rec.type==2] = cpn.m[rec.type==2]
 for (i in which(rec.type==3)) {
  rec.bp = c(ceiling(runif(1,0,1)*cpn.p[i]),ceiling(runif(1,0,1)*cpn.m[i]))
  if (rec.bp[2] < cpn.m[i]) {temp = c(pop.p[i,1:rec.bp[1]],pop.m[i,(1+rec.bp[2]):cpn.m[i]])}
  else {temp = pop.p[i,1:rec.bp[1]]}
  
  if (length(temp) <= ncol(pop.o)) {pop.o[i,1:length(temp)] = temp}
  else {
   pop.o = cbind(pop.o,matrix(0,2*ne.n,(length(temp)-ncol(pop.o))))
   pop.o[i,] = temp
  }
  cpn.o[i] = length(temp)
 }
 for (i in which(rec.type==4)) {
  rec.bp = c(ceiling(runif(1,0,1)*cpn.m[i]),ceiling(runif(1,0,1)*cpn.p[i]))
  if (rec.bp[2] < cpn.p[i]) {temp = c(pop.m[i,1:rec.bp[1]],pop.p[i,(1+rec.bp[2]):cpn.p[i]])}
  else {temp = pop.m[i,1:rec.bp[1]]}
  
  if (length(temp) <= ncol(pop.o)) {pop.o[i,1:length(temp)] = temp}
  else {
   pop.o = cbind(pop.o,matrix(0,2*ne.n,(length(temp)-ncol(pop.o))))
   pop.o[i,] = temp
  }
  cpn.o[i] = length(temp)
 }
 
 nmut = rpois(1,2*ne.n*ncol(pop.o)*mu)
 if (nmut > 0) {
  for (i in 1:nmut) {
   lmut = c(sample(1:(2*ne.n),1),sample(1:ncol(pop.o),1))
   if (pop.o[lmut[1],lmut[2]] > 0) {
    allele = ceiling(runif(1,0,10000))
    pop.o[lmut[1],lmut[2]] = allele
   }
  }
 }
 
 x = rep(0,ne.n)
 y = rep(0,ne.n)
 fit = rep(0,ne.n)
 for (i in 1:ne.n) {
  temp = c(pop.o[i,],pop.o[ne.n+i,])
  y[i] = sum(temp>0)
  x[i] = sum(unique(temp)>0)
  if (y[i] > ncol(fitlist)) {fit[i] = 0}
  else {fit[i] = fitlist[x[i],y[i]]}
 }
  
 return(list(copies = pop.o, cpn = cpn.o, fit = fit, x = x, y = y))

}

### Function to migrate individuals from one population to another (only one direction; needed to be run twice if migration is reciprocal

mhc.mig = function (pop1, pop2, m12) {#pop1 = donor population, pop2 = recipient population, m12 = probability of an individual from pop1 will migrate

 ne.1 = length(pop1$cpn)/2
 ne.2 = length(pop2$cpn)/2
 coldiff = ncol(pop2$copies) - ncol(pop1$copies)
 pop2.o = pop2
 if (coldiff < 0) {pop2.o$copies = cbind(pop2.o$copies,matrix(0,2*ne.2,abs(coldiff)))}

 temp = which(rbinom(ne.1,1,m12)==1)
 if (length(temp) > 0) {
  temp.mp = pop1$copies[temp,]
  temp.mm = pop1$copies[(ne.1+temp),]
  if (length(temp) == 1) {
   temp.mp = t(as.matrix(temp.mp))
   temp.mm = t(as.matrix(temp.mm))
  }
  if (coldiff > 0) {
   temp.mp = cbind(temp.mp, matrix(0,length(temp),coldiff))
   temp.mm = cbind(temp.mm, matrix(0,length(temp),coldiff))
  }
  
  pop2.o$copies = rbind(pop2.o$copies[1:ne.2,],temp.mp,pop2.o$copies[(ne.2+1):(2*ne.2),],temp.mm)
  pop2.o$cpn = c(pop2.o$cpn[1:ne.2],pop1$cpn[temp],pop2.o$cpn[(ne.2+1):(2*ne.2)],pop1$cpn[ne.1+temp])
  pop2.o$fit = c(pop2.o$fit[1:ne.2],pop1$fit[temp])
  pop2.o$x = c(pop2.o$x[1:ne.2],pop1$x[temp])
  pop2.o$y = c(pop2.o$y[1:ne.2],pop1$y[temp])
 
 }
 return (pop2.o) #Only the new pop2 is returned

}

### Function to calculate statistics of a population and put it into a vector for output

mhc.desc = function (pop) { #Describe a population
 #The statistics in order: Ne Mean_x SD_x Min_x Max_x Mean_y SD_y Min_y Max_y Mean_ratio SD_ratio Min_ratio Max_ratio Mean_fit SD_fit Min_fit Max_fit Mean_hapuniq SD_hapuniq Min_hapuniq Max_hapuniq Num_alleles Effnum_alleles Freq_1st Freq_2nd Freq_3rd Freq_4th Freq_5th

 out = NULL
 out[1] = nrow(pop$copies)/2
 out[2] = mean(pop$x)
 out[3] = sd(pop$x)
 out[4:5] = range(pop$x)
 out[6] = mean(pop$y)
 out[7] = sd(pop$y)
 out[8:9] = range(pop$y)
 out[10] = mean(pop$x/pop$y)
 out[11] = sd(pop$x/pop$y)
 out[12:13] = range(pop$x/pop$y)
 out[14] = mean(pop$fit)
 out[15] = sd(pop$fit)
 out[16:17] = range(pop$fit)
 haplo = apply(pop$copies,1,uniqnum)
 out[18] = mean(haplo)
 out[19] = sd(haplo)
 out[20:21] = range(haplo)
 allst = as.numeric(table(as.numeric(pop$copies)))[-1]
 out[22] = length(allst)
 out[23] = 1/sum((allst/sum(allst))^2)
 
 temp = sort(table(as.numeric(pop$copies)),decreasing=TRUE)[2:6]
 temp2 = as.numeric(names(temp))
 out[24] = sum(apply(pop$copies==temp2[1],1,sum)>0)
 out[25] = sum(apply(pop$copies==temp2[2],1,sum)>0)
 out[26] = sum(apply(pop$copies==temp2[3],1,sum)>0)
 out[27] = sum(apply(pop$copies==temp2[4],1,sum)>0)
 out[28] = sum(apply(pop$copies==temp2[5],1,sum)>0)
 return(out)
}


### Commands needed for variable recombination rate

pop$rrm = rep(1,2*ne) #Initialize rrm as 1 for all chromosomes
 
 #Run every generation: mutate RRM with probability 0.002 to increase by exp(0.05) and same to decrease
{
 temprr = sample(c(exp(-0.05),1,exp(0.05)),2*ne,replace=TRUE,prob=c(0.002,0.996,0.002))
 pop$rrm = pop$rrm*temprr
}

### Alternate function to evolve a population, with recombination rate variable

mhc.evolve.rrate = function (pop, fitlist, rec, mu, ne.n) {#The rec rate in the input is "baseline"; the chromosomal recombination rate is recorded inside the population data as pop$rrm
 ne = length(pop$cpn)/2
 par1 = sample(1:ne,ne.n,replace=TRUE,prob=pop$fit) #Choose parent individuals (diploid)
 par2 = sample(1:ne,ne.n,replace=TRUE,prob=pop$fit) 
 
 rec.true = rec*pop$rrm #Recombination rate by chromosome, equals baseline * multiplier
 rec.true2 = sqrt(rec.true[1:ne] * rec.true[(ne+1):(2*ne)]) #Recombination rate by individual, geometric mean of both chromosomes
 rec.type.yn = rbinom(2*ne.n, 1, rec.true2[c(par1,par2)]) #Whether each offspring chromosome is a recombinant
 rec.type = rec.type.yn*2 + sample(1:2,2*ne.n,replace=TRUE) #0*2+1 = 1, 0*2+2 = 2, 1*2+1 = 3, 1*2+2 = 4. The rest shall be the same.
 rec.new.p = pop$rrm[c(par1,par2)] #RRM of the offspring generation, inherited from the "head" chromosome.
 rec.new.m = pop$rrm[c(par1+ne,par2+ne)]
 rec.new = rep(0,2*ne.n)
 rec.new[which(rec.type%%2==1)] = rec.new.p[which(rec.type%%2==1)]
 rec.new[which(rec.type%%2==0)] = rec.new.m[which(rec.type%%2==0)]
 
 #rec.help = c((1-rec)/2,(1-rec)/2,rec/2,rec/2)
 #rec.type = sample(1:4,2*ne.n,replace=TRUE,prob=rec.help) # 1-> Gamete from paternal copy, 2-> Gamete from maternal copy, 3-> Paternal first half, maternal second half, 4-> vv
 
 pop.p = pop$copies[c(par1,par2),] #Paternal copy in father, then paternal copy in mother (for meiotic recombination)
 pop.m = pop$copies[c(par1+ne,par2+ne),]
 cpn.p = pop$cpn[c(par1,par2)]
 cpn.m = pop$cpn[c(par1+ne,par2+ne)]
 
 pop.o = matrix(0,2*ne.n,ncol(pop$copies))
 cpn.o = rep(0,2*ne.n)
 pop.o[rec.type==1,] = pop.p[rec.type==1,]
 pop.o[rec.type==2,] = pop.m[rec.type==2,]
 cpn.o[rec.type==1] = cpn.p[rec.type==1]
 cpn.o[rec.type==2] = cpn.m[rec.type==2]
 for (i in which(rec.type==3)) {
  rec.bp = c(ceiling(runif(1,0,1)*cpn.p[i]),ceiling(runif(1,0,1)*cpn.m[i]))
  if (rec.bp[2] < cpn.m[i]) {temp = c(pop.p[i,1:rec.bp[1]],pop.m[i,(1+rec.bp[2]):cpn.m[i]])}
  else {temp = pop.p[i,1:rec.bp[1]]}
  
  if (length(temp) <= ncol(pop.o)) {pop.o[i,1:length(temp)] = temp}
  else {
   pop.o = cbind(pop.o,matrix(0,2*ne.n,(length(temp)-ncol(pop.o))))
   pop.o[i,] = temp
  }
  cpn.o[i] = length(temp)
 }
 for (i in which(rec.type==4)) {
  rec.bp = c(ceiling(runif(1,0,1)*cpn.m[i]),ceiling(runif(1,0,1)*cpn.p[i]))
  if (rec.bp[2] < cpn.p[i]) {temp = c(pop.m[i,1:rec.bp[1]],pop.p[i,(1+rec.bp[2]):cpn.p[i]])}
  else {temp = pop.m[i,1:rec.bp[1]]}
  
  if (length(temp) <= ncol(pop.o)) {pop.o[i,1:length(temp)] = temp}
  else {
   pop.o = cbind(pop.o,matrix(0,2*ne.n,(length(temp)-ncol(pop.o))))
   pop.o[i,] = temp
  }
  cpn.o[i] = length(temp)
 }
 
 nmut = rpois(1,2*ne.n*ncol(pop.o)*mu)
 if (nmut > 0) {
  for (i in 1:nmut) {
   lmut = c(sample(1:(2*ne.n),1),sample(1:ncol(pop.o),1))
   if (pop.o[lmut[1],lmut[2]] > 0) {
    allele = ceiling(runif(1,0,10000))
    pop.o[lmut[1],lmut[2]] = allele
   }
  }
 }
 
 x = rep(0,ne.n)
 y = rep(0,ne.n)
 fit = rep(0,ne.n)
 for (i in 1:ne.n) {
  temp = c(pop.o[i,],pop.o[ne.n+i,])
  y[i] = sum(temp>0)
  x[i] = sum(unique(temp)>0)
  if (y[i] > ncol(fitlist)) {fit[i] = 0}
  else {fit[i] = fitlist[x[i],y[i]]}
 }
  
 return(list(copies = pop.o, cpn = cpn.o, fit = fit, x = x, y = y, rrm = rec.new))

}


### Example 1: burn-in for 20,000 generations

ne = 1000 #Population size is 1000 diploid individuals
sel_str = 2 #Intermediate selection level, with s_x = 0.02, s_y = 0.005
outfile = "Output_report.txt"
outfile2 = "Output.txt"
#Placeholder file names

s1 = 0.01*sel_str
s2 = 0.0025*sel_str
k1 = 0.95
k2 = 1.05
fitlist = matrix(0,100,100)
fitlist[1,2] = (1+s1)
fitlist[2,2] = (1+s1)^(1+k1)
for (x in 1:100) {
 for (y in 3:100) {
  if (x <= y) {
   t1 = (1+s1)^(sum(k1^(0:(x-1))))
   t2 = (1-s2)^(sum(k2^(0:(y-3))))
   fitlist[x,y] = t1*t2
  }
 }
}

pop = list(copies = matrix(1,2*ne,5), cpn = rep(5,2*ne), fit = rep(1,ne), x = rep(1,ne), y = rep(10,ne))

for (gen in 1:20000) {
 pop = mhc.evolve(pop,fitlist,0.01,0.0005,ne)
 if (gen %% 100 == 0) {
  write(c(mean(pop$x),sd(pop$x),mean(pop$y),sd(pop$y),mean(pop$fit),sd(pop$fit)),outfile,ncolumns=6,append=TRUE)
  #This writes to the "report" file each 100 generations, so that one can verify that values such as population mean x, y and fitness have reached equilibrium.
 }

}

write.table(pop$copies,outfile2,col.names=FALSE,row.names=FALSE)
#This writes to the population snapshot/data file. In this case, the 1st-1000th rows are the "paternal" chromosomes of each individual, and 1001st-2000th rows are "maternal" chromosomes. Within each row the alleles are ordered and "extra empty space" marked with 0s.


### Example 2: reciprocal migration

ne = 1000 #Two populations are involved and the size of EACH is 1000
sel_str = 2 #Intermediate selection level, with s_x = 0.02, s_y = 0.005
mig = 1 #Migration rate is 1 individual (on average) per generation per direction

infile1 = "Input_pop1.txt"
infile2 = "Input_pop2.txt"
outfile1 = "Output_pop1.txt"
outfile2 = "Output_pop2.txt"
outana = "Output_report.txt"
#Placeholder file names

s1 = 0.01*sel_str
s2 = 0.0025*sel_str
k1 = 0.95
k2 = 1.05
fitlist = matrix(0,100,100)
fitlist[1,2] = (1+s1)
fitlist[2,2] = (1+s1)^(1+k1)
for (x in 1:100) {
 for (y in 3:100) {
  if (x <= y) {
   t1 = (1+s1)^(sum(k1^(0:(x-1))))
   t2 = (1-s2)^(sum(k2^(0:(y-3))))
   fitlist[x,y] = t1*t2
  }
 }
}

pop1 = mhc.read(infile1,fitlist)
pop2 = mhc.read(infile2,fitlist)

txt = "Gen Pop Ne Mean_x SD_x Min_x Max_x Mean_y SD_y Min_y Max_y Mean_ratio SD_ratio Min_ratio Max_ratio Mean_fit SD_fit Min_fit Max_fit Mean_hapuniq SD_hapuniq Min_hapuniq Max_hapuniq Num_alleles Effnum_alleles Freq_1st Freq_2nd Freq_3rd Freq_4th Freq_5th"
write(txt,outana)
#This is the header for the report file, which would be updated each 10 generations. The report file can be read in R with header=TRUE.

for (gen in 1:10000) {
 pop1x = mhc.mig(pop2, pop1, mig/ne)
 pop2x = mhc.mig(pop1, pop2, mig/ne)
 #These temporary populations are derived from reciprocal migration
 
 pop1 = mhc.evolve(pop1x,fitlist,0.01,0.0005,ne)
 pop2 = mhc.evolve(pop2x,fitlist,0.01,0.0005,ne)
 #And then evolved to produce the next generation
 
 if (gen %% 10 == 0) {
  out1 = mhc.desc(pop1)
  write(c(gen,1,out1),outana,ncolumns=30,append=TRUE)
  out2 = mhc.desc(pop2)
  write(c(gen,2,out2),outana,ncolumns=30,append=TRUE)
  #Each 10 generations the population statistics are written to the report file, with pop1 and pop2 alternating.
 }
 

}

write.table(pop1$copies,outfile1,col.names=FALSE,row.names=FALSE)
write.table(pop2$copies,outfile2,col.names=FALSE,row.names=FALSE)
#Both populations at end of run also written to files.
