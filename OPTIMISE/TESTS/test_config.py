#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 14:45:58 2023

@author: andrea bonato
"""
import sys
import os
path = os.getcwd()
sys.path.insert(0, path+'/..')
import functions


success = functions.read_config("opti_input_example.txt")