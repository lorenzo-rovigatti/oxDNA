#Code to map a multivariate normal distribution with initial mean and covariance
#to a target multivariate normal distribution.
#It minimises the relative entropy of the two distributions by using gradient descent techniques 
#integrated with reweighting of the mean and covariance

import numpy as np
import math
import copy
from scipy import optimize

dimension = 2
par_dimension = 6
#tolerance = 0.001


#Mean and covariance of the target distribution
target_cov = np.zeros([dimension,dimension], dtype=float)
target_mu = np.zeros(dimension, dtype=float)

target_cov[0,0] = 2.5
target_cov[0,1] = 0.6
target_cov[1,0] = 0.6
target_cov[1,1] = 1.3
target_mu[0] = 1.
target_mu[1] = 1.3

#print(target_cov.trace())
#print(np.linalg.det(target_cov))

target_inv_cov = np.linalg.inv(target_cov)
det_target_inv_cov = np.linalg.det(target_inv_cov)

data = [] #stores sampled data

#maps the covaraince matrix and mean to a vector
def cov_and_mean_to_vector(cov,mu):
    
    par_v = []
    
    for i in range(dimension) :
        for j in range(dimension) :
            par_v.append(cov[i][j])
        
    for i in range(dimension) :
        par_v.append(mu[i])
        
    return par_v

#maps a vector to the covaraince matrix and mean
def vector_to_cov_and_mean(par_v):
    
    cov = np.zeros([dimension,dimension], dtype=float)    
    mu = np.zeros(dimension, dtype=float)
    
    for i in range(dimension) :
        for j in range(dimension) :
            cov[i][j] = par_v[i*dimension+j]
        
    for i in range(dimension) :
        mu[i] = par_v[dimension*dimension+i]
        
    return cov,mu


#maps the covariance matrix, the mean and a parameter vector to a vector
#this will be used for the oxDNA version of the optimisation procedure
def encode_args(cov,mu,par):
    
    args = np.zeros(dimension*(dimension+1)+par_dimension,dtype=float)
    
    for i in range(dimension) :
        for j in range(dimension) :
            args[i*dimension+j] = cov[i][j]
        
    for i in range(dimension) :
        args[dimension*dimension+i] = mu[i]
        
    for i in range(par_dimension) :
        args[dimension*(dimension+1)+i] = par[i]
        
    return args


#inverse of encode_argv
def decode_args(args):
    
    cov = np.zeros([dimension,dimension], dtype=float)    
    mu = np.zeros(dimension, dtype=float)
    par = np.zeros(par_dimension,dtype=float)
    
    for i in range(dimension) :
        for j in range(dimension) :
            cov[i][j] = args[i*dimension+j]
        
    for i in range(dimension) :
        mu[i] = args[dimension*dimension+i]
        
    for i in range(par_dimension) :
        par[i] = args[dimension*(dimension+1)+i]
        
    return cov,mu,par

#Sample from a multivariate normal distribution with given mean and covariance (stored in par, see cov_and_mean_to_vector)
#if save = True, the sampled data is stored inside the global list data.
#returns mean and covariance of the sample
def normal_sample_average(par,save=True) :
    
    cov,mu = vector_to_cov_and_mean(par)
    #print(mu)
    #print(cov)
    global data
    
    data_al = np.random.multivariate_normal(mu,cov,100000).T
    
    if save :
        data = data_al    
    
    mean = np.zeros(dimension,dtype=float)
    covariance = np.zeros((dimension,dimension), dtype=float)
    
    #compute mean
    for i in range(len(data_al[0])) :
        mean[0]+=data_al[0,i]/(1.0*len(data_al[0]))
        mean[1]+=data_al[1,i]/(1.0*len(data_al[0]))
        
    #compute covariance
    for i in range(len(data_al[0])) :
        for j in range(dimension) :
            for k in range(j,dimension) :
                covariance[j,k]+=(data_al[j][i]-mean[j])*(data_al[k][i]-mean[k])/(1.0*len(data_al[0]))
                
    for j in range(dimension) :
        for k in range(j+1,dimension) :
                covariance[k,j] = covariance[j,k]
    
    return covariance,mean

#Compute mean and covariance with parameters par, by reweighting mean and covariance with parameters par0.
#The data used to compute ensamble avarages is stored inside the global list data.
#Data must be sampled with parameters par0 before reweighting.
#The reweighting procedure is for multivariate normal distributions.
def normal_reweight_cov_and_mu(par,par0) :
    
    cov0,mu0 = vector_to_cov_and_mean(par0)    
    cov1,mu1 = vector_to_cov_and_mean(par)
    
    cov = np.zeros([dimension,dimension], dtype=float)    
    mu = np.zeros(dimension, dtype=float)
    
    #sample = np.random.multivariate_normal(mu,cov,100000).T
    
    #<e^-DH>
    av_e_to_deltaH = 0.

    #data are sampled with par0 before calling reweighting, and stored globally.
    #This is to avoid sampling multiple times when unnecessary.
    #Fits well with what we want to do with oxDNA
    
    #reweight mean
    for i in range(len(data[0])) :
        
        x = np.array([data[0,i],data[1,i]])
        
        v0 = x-mu0
        v1 = x-mu1
        
        #This is for a multivariate normal distribution
        deltaH = +0.5*np.dot(v1.transpose(),np.dot(np.linalg.inv(cov1),v1)) +\
            -0.5*np.dot(v0.transpose(),np.dot(np.linalg.inv(cov0),v0))
        
        mu[0]+=data[0,i]*math.exp(-deltaH)
        mu[1]+=data[1,i]*math.exp(-deltaH)
        
        av_e_to_deltaH += math.exp(-deltaH)
        
    mu = mu*(1./av_e_to_deltaH)
    
    #reweight covariance
    for i in range(len(data[0])) :
        for j in range(dimension) :
            for k in range(j,dimension) :
                
                x = np.array([data[0,i],data[1,i]])
                
                v0 = x-mu0
                v1 = x-mu1
                
                deltaH = 0.5*np.dot(v1.transpose(),np.dot(np.linalg.inv(cov1),v1)) +\
                    -0.5*np.dot(v0.transpose(),np.dot(np.linalg.inv(cov0),v0))
                
                cov[j,k]+=(data[j][i]-mu[j])*(data[k][i]-mu[k])*math.exp(-deltaH)/(1.0*av_e_to_deltaH)
                
    for j in range(dimension) :
        for k in range(j+1,dimension) :
                cov[k,j] = cov[j,k]
    
    return cov,mu


#Compute Relative Entropy.
#It is a function of some parameters par;
#args are extra parameters:
#the mean and covariance of a given set of parameters par0, and par0 itself.
#If par = par0, the function uses the mean and covariance in args to compute the relative entropy,
#otherwise it estimates mean and covariance for par by reweighting par0.
#This function is built this way so that optimize.minimize uses reweighting to estimate Rel_entropy(par0+delta_par)
#e.g. when estimating the gradient and Hessian of the Relative entropy in par0
def Relative_entropy_wRew(par,args):
    
    cov0,mu0,par0 = decode_args(args)
    rew = False
    
    #check if par == par0 (1/1000 tolerance)
    for i in range(par_dimension) :
        if par[i] < par0[i]*(999./1000.) or par[i] > par0[i]*(1001./1000.) :
            rew = True
            break
    
    #if par != par0, then reweight mean and covariance from mean and covariance with par0
    if rew :
        cov,mu = normal_reweight_cov_and_mu(par,par0)
    else :
        cov,mu = cov0,mu0
    

    #print(cov)
    #print(mu)
    
    #compute relative entropy
    
    S = 0.
    
    S += np.dot(cov,target_inv_cov).trace()
    
    S -= math.log(det_target_inv_cov*np.linalg.det(cov))+dimension
    
    delta_mu = mu - target_mu
    
    S += np.dot(np.dot(delta_mu.transpose(),target_inv_cov),delta_mu)
    
    S = S*0.5
    
    return S


par0 = [2.8,0.5,0.5,1.0,0.7,1.8]


"""
#test reweighting

par1 = [3.3,0.57,0.57,2.2,0.6,1.7]

cov,mu = normal_sample_average(par0,True)

print(cov)
print(mu)

cov,mu = normal_sample_average(par1,False)

print(cov)
print(mu)

cov,mu = normal_reweight_cov_and_mu(par1,par0)

print(cov)
print(mu)
"""

#optimise

#1 step sampling, 3 steps pure reweighting



nsteps = 0 #number of steps

par = par0

while nsteps < 60 : #impose tolerance here!

    #sample covariance and mean
    cov,mu = normal_sample_average(par,True)
    arguments = encode_args(cov,mu,par)
    
    sol = optimize.minimize(Relative_entropy_wRew,par,args=(arguments),method='nelder-mead',options={'maxiter':3})
    
    par = sol.x
    
    #Ad hoc: symmetrize covariance (not needed for oxDNA3)
    #right procedure: use minimize with constraints!!    
    av_diag = (par[1]+par[2])/2.
    
    par[1] = av_diag
    par[2] = av_diag
    
    print(sol)
    
    nsteps += 3
    
    print(str(nsteps)+": "+str(par))