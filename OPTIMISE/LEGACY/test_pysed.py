#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 12:42:58 2024

@author: yqb22156

"""

import pysed



file = "./tets_pysed.txt"
strold = "T = *K"    
strnew = "T = "+str(T)+"K"    

pysed.replace(strold,strnew,file)   