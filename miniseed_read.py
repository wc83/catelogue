#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 10:30:07 2018

@author: william
"""

from obspy import read
stre = read("/Users/william/Documents/scanner/output_data/EXP_data_stream_v1.mseed")
print(stre)
