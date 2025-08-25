#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 11:59:33 2024

@author: yqb22156
"""

import pysed

seq_id = 2
file_name = "input_VMMC"

def set_topo_input_file(seq_id,file_name) :
    
    strold1 = r"(topology) = .*(top)$"
    strnew = "topology = topo_s"+str(seq_id)+".top"
    
    pysed.replace(strold1,strnew,file_name) 
    
    return

set_topo_input_file(seq_id,file_name)