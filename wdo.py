#!/usr/bin/python
""" --------------------------------------------------------------
Sample Python Code for the Wind Driven Optimization.
Optimization of the Sphere Function in the range of [-5, 5].
by Zikri Bayraktar, PhD. - thewdoalgorithm@gmail.com

DISCLAIMER: This code is provided for educational purposes only. 
USE IT AT YOUR OWN RISK! Code is by no means a highly optimized code,
hence do not use it to compare other algorithms without studying the
algorithm and tuning for your purposes.

This script is developed in Anaconda Spyder v.2.3.4 running Python 2.7.9
-------------------------------------------------------------------------
Please refer to the following journal article in your research papers:
Z. Bayraktar, M. Komurcu, J. A. Bossard and D. H. Werner, "The Wind 
Driven Optimization Technique and its Application in Electromagnetics," 
IEEE Transactions on Antennas and Propagation, Volume 61, Issue 5, 
pages 2745 - 2757, May 2013.
--------------------------------------------------------------- """

# import some modules:
import math
import random
from datetime import datetime
import numpy as np # import NumPy for numeric calculations.
import scipy as sp 
import matplotlib.pyplot as plt

random.seed(7)  # initialize a random seed for repeatability of pseudorandom numbers

# Retrieve time/date for log purposes:
now = datetime.now()
current_year = now.year
current_month = now.month
current_day = now.day
current_hour = now.hour
current_minute = now.minute
current_second = now.second

# Print variables named current_*** for log purposes: 
print '%s-%s-%s' % (current_year, current_month, current_day)
print '%s:%s:%s' % (current_hour, current_minute, current_second)

# WDO related variable settings:
popsize = 20       	# population size
npar = 5			# dimension of the problem
maxit = 500			# max number of iterations
RT = 3				# RT coefficient
g = 0.20			# gravitational constant
alp = 0.4			# constants in the update eq.
c = 0.4				# coriolis effect
maxV = 0.3			# max allowed speed
dimMin = -5			# lower dimension boundary
dimMax = 5			# upper dimension boundary
xmin = 0.0			# minimum x value per dimension
#----------------------------------------------------------

# Initialize WDO population position, and velocity:
# Randomize population in the range of [-1,1]:
pos = 2 * (np.random.rand(popsize, npar) - 0.5)
vel = maxV * 2 * (np.random.rand(popsize, npar) - 0.5)
#----------------------------------------------------------

# Evaluate initial population: (Sphere Function)
pres = np.zeros((popsize))  # initialize the pres vector to zeros.
a=pres.shape # size of the pres

for i in range(0, popsize): #pyhon index starts from 0 !!!
    x = ((dimMax-dimMin) * (pos[i,:]+1)/2) + dimMin
    temp = np.power((x-xmin),2)
    pres[i] = temp.sum()

#----------------------------------------------------------
# Finding best air parcel in the initial population :

globalpres = pres.min()         # minimum pressure
minIndx = np.where(pres == pres.min())  # index of minimum pressure
globalpos = pos[minIndx,:]       # position vector for the minimum

minpres = np.zeros(maxit)
keepglob = np.zeros(maxit)

indx = np.argsort(pres)      # index of sorted
pos = pos[indx,:]
minpres[0] = globalpres		# save the minimum pressure
keepglob[0] = globalpres;   
#-----------------------------------------------------------------

velot = np.zeros((popsize, npar))
keepglob = np.ones(maxit)
	
# Start iterations:
itr = 1 #iteration counter
for ij in range(1,maxit):
    #update velocity
    for i in range(popsize):
        a = np.random.permutation(range(0,npar)) #random perm    
        velot[i,:] = 1*vel[i,a]
        vel[i,:] =  (1-alp) * vel[i,:] - (g*pos[i,:]) + \
                    abs((1/(i+1))-1) *((globalpos-pos[i,:]*RT)) + \
                    (c*velot[i,:]/(i+1))
        #python index starts from zero, watch out for division by zero error!

    #check velocity
    vel=vel.clip(-maxV,maxV)
    #update air parcel position
    pos=pos+vel
    pos=pos.clip(-1,1)
    
    #evaluate the new position
    for i in range(0, popsize): #pyhon index starts from 0 !!!
        x = ((dimMax-dimMin) * (pos[i,:]+1)/2) + dimMin
        temp = np.power((x-xmin),2)
        pres[i] = temp.sum()

    # Finding best air parcel in the initial population :

    mpres = pres.min()         # minimum pressure
    mIndx = np.where(pres == pres.min())  # index of minimum pressure
    gpos = pos[mIndx,:]       # position vector for the minimum

    indx = np.argsort(pres)     # index of sorted
    pos = pos[indx,:]           #sort position
    vel = vel[indx,:]           #sort velocity
    
    if mpres < globalpres:  #if lower pressure found, update the global min
        globalpres = mpres
        globalpos = gpos
        
    keepglob[ij] = globalpres
    
    #minpres[ij] = globalpres		# save the minimum pressure
    #keepglob[ij] = globalpres;   

    
plt.plot(keepglob)
plt.ylabel('Pressure')
plt.xlabel(('Iteration'))
plt.yscale('log', nonposy='clip')
plt.show()
    
#----------------------------------------------------------
