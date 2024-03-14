#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:19:53 2024

@author: yqb22156
"""

from scipy import optimize
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


a = []
for i in range(size):
    a.append(1.)


def parallel_function_caller(x,stopp):
    stopp[0]=comm.bcast(stopp[0], root=0)
    summ=0
    if stopp[0]==0:
        #your function here in parallel
        x=comm.bcast(x, root=0)
        loc_term = (x[0]-rank-1)*(x[0]-rank-1)*a[rank]
        summ=comm.reduce(loc_term,op=MPI.SUM, root=0)
        if rank==0:
            print(x[0],summ)
    return summ

if rank == 0 :
   stop=[0]
   x = np.zeros(1)
   x[0]=20
   #xs = minimize(parallel_function_caller, x, args=(stop))
   xs = optimize.minimize(parallel_function_caller,x0= x, args=(stop,))
   print ("the argmin is "+str(xs))
   stop=[1]
   parallel_function_caller(x,stop)

else :
   stop=[0]
   x=np.zeros(1)
   ite = 0
   while stop[0]==0:
      ite += 1
      summ = parallel_function_caller(x,stop)
      print("while: rank "+ str(rank)+ ", ite "+ str(ite)+ ", x " + str(x[0]) + ", summ " + str(summ))