#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 16:06:30 2018

@author: william
"""

def find_nearest(array,value):
    import numpy as np
    idx = (np.abs(array-value)).argmin()
    return (array[idx],idx)