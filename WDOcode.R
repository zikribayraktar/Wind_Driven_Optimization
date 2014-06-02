#-------------------------------------------------
# Sample R Code for the Wind Driven Optimization.
# Optimization of the Sphere Function in the range of [-5, 5].
# by Dr. Zikri Bayraktar - thewdoalgorithm@gmail.com
#
# DISCLAIMER: This code is provided for educational purposes
# only. Use it at your own risk!
#-------------------------------------------------
#
# Please refer to the following journal article in your publications:
# Z. Bayraktar, M. Komurcu, J. A. Bossard and D. H. Werner, "The Wind 
# Driven Optimization Technique and its Application in Electromagnetics," 
# IEEE Transactions on Antennas and Propagation, Volume 61, Issue 5, 
# pages 2745 - 2757, May 2013.
#-------------------------------------------------
#
rm(list=ls(all=TRUE))             #clear environment from memory
cat("\014")                       #clear console
starttime <- Sys.time()           #record start time
(WD <- getwd())                   #findout the current directory
if (!is.null(WD)) setwd(WD)       #set to working directory
if(file.exists("WDOout.txt")) file.remove("WDOout.txt") #remove old files
if(file.exists("WDOprs.txt")) file.remove("WDOprs.txt")
if(file.exists("WDOpos.txt")) file.remove("WDOpos.txt")
#---------------------------------
file.create("WDOout.txt")      #create new files
file.create("WDOprs.txt")
file.create("WDOpos.txt")
#---------------------------------
# User defined WDO parameters
set.seed(007)                     # set the seed for random # generator
param.popsize <- 9
param.npar <- 4
param.maxit <- 500
param.RT <- 3
param.g <- 0.2
param.alp <- 0.4
param.c <- 0.4
param.maxV <- 0.3
param.dimMin <- -5
param.dimMax <- 5
pres <- vector()      # to record the pressure per member
minpres <- vector()   # to record the min pressure per iteration
minpos <- matrix(data=NA,nrow=param.maxit,ncol=param.npar)    # to record the position with min pressure per iteration
velot <- vector()
globpos <- vector()
velnew <- vector()
keeppres <- vector()
keepindx <- vector()
#---------------------------------
# Initialize WDO population, position and velocity:
# Randomize population in the range of [-1, 1]:
pos <- matrix( (2 * (runif(param.popsize * param.npar) -0.5)), ncol=param.npar) 
# Randomize the velocity in the range of [-maxV maxV]
vel <- matrix( (param.maxV * 2 * (runif(param.popsize * param.npar) -0.5)), ncol=param.npar) 

# Evaluate initial population: (Sphere Function)
for ( i in 1:param.popsize)
{
  x <- (param.dimMax - param.dimMin) * ((pos[i,1:param.npar]+1)/2) + param.dimMin;
  pres[i] <- sum(x*x)
  #print(pres[i])
}

# Sort and find the minimum and its index for the initial evaluation:
sorted <- sort(pres, index.return=TRUE)
newpos <- pos[sorted$ix, 1:param.npar] 
minpres[1] <- min(sorted$x)
globpres <- minpres[1]
minpos[1,1:param.npar] <- newpos[1,1:param.npar]
pos <- newpos
globpos <- minpos[1,1:param.npar]

# Generate an index from 1 to param.popsize. This will be used to select other vel dim.
xs <- 1:param.npar

# Start the iterations:
iter <- 1
for ( ij in (2:param.maxit) )
{

  #update the velocity of each member:
  for ( jj in (1:param.popsize) ) 
  {
  # randomize dimensions for velocity update:
    ys <- sample(xs)
    velot <- vel[jj,ys]    
    vel[jj,(1:param.npar)] <- ((1-param.alp) * vel[jj,(1:param.npar)]) - ((param.g*pos[jj,(1:param.npar)])) + abs(1-(1/jj))*((globpos-pos[jj,(1:param.npar)])*param.RT) + param.c*(velot/jj)
  }
  
  # Check and limit velocity:
  vel[(1:param.popsize), (1:param.npar)] <- pmax(vel[(1:param.popsize), (1:param.npar)], -param.maxV)
  vel[(1:param.popsize), (1:param.npar)] <- pmin(vel[(1:param.popsize), (1:param.npar)], param.maxV)

  # update air parcel position then limit:
  pos <- pos + vel
  pos[(1:param.popsize), (1:param.npar)] <- pmax(pos[(1:param.popsize), (1:param.npar)], -1)
  pos[(1:param.popsize), (1:param.npar)] <- pmin(pos[(1:param.popsize), (1:param.npar)], 1)

  # Evaluate population: (Sphere Function)
  for ( i in (1:param.popsize))
  {
    x <- (param.dimMax - param.dimMin) * ((pos[i,1:param.npar]+1)/2) + param.dimMin;
    pres[i] <- sum(x*x)
  }
  
  
  # Sort and find the minimum and its index for this iterations:
  sorted <- sort(pres, index.return=TRUE)
  newpos <- pos[sorted$ix, (1:param.npar)] 
  newvel <- vel[sorted$ix, (1:param.npar)] 
  minpres[ij] <- min(sorted$x)

  minpos[ij,(1:param.npar)] <- newpos[1,(1:param.npar)]
  
  pos <- newpos
  vel <- newvel
  
  if (minpres[ij] < globpres)
  {  globpres <- minpres[ij]
     globpos <- minpos[ij,1:param.npar]
     print(c(ij, globpres))
  }
  keeppres[ij] <- globpres
  keepindx[ij] <- ij
  
} #(ij) - end-of-iterations

  plot(keepindx, keeppres, type="l", xlab="Iteration",ylab="Min Pressure", log="y")

#---------------------------------
endtime <- Sys.time()
runtime <- endtime - starttime
save(list=ls(all=TRUE),file="SaveSession")
#to recover your session:
#load("SaveSession")
# end-of-file