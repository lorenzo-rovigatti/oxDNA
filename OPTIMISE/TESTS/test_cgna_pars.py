#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 11:16:58 2023

@author: andrea bonato
"""
import sys
import os
path = os.getcwd()
sys.path.insert(0, path+'/..')
import get_cgdna_pars as gcp
import numpy as np


#test

ids = [1,6,7,8,11]

N = 20

bp_ave_gs = np.load(gcp.path_averages+"cgna_ave_gs.npy")
bp_ave_cov = np.load(gcp.path_averages+"cgna_ave_cov.npy")

gs, cov = gcp.get_target_mean_and_covariance_diag_ave(N, ids)

 