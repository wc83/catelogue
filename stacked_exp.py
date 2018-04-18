#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 09:45:52 2018

@author: william
"""

def stacked_exp(fmin.fmax):
    import numpy as np
    from obspy.core import read
    #%%
    
    
    
    seisb = read('/Users/william/Documents/lb01/14_326z.mseed')
    stream=seisb.copy()
    stream=seisb.clear()
    stream2=stream.copy()
    #%%
    for x in range(0,10,1):
        if x==0:
            seis = read('/Users/william/Documents/lb01/14_326z.mseed')
            M=[[1,36,46],[2,37,17],[4,44,34],[12,19,2],[20,21,27]]
        if x==1:
            seis = read('/Users/william/Documents/lb01/14_327z.mseed')
            M=[[3,3,40],[5,59,54],[9,47,32],[15,18,34],[18,39,4]]
        if x==2:
            seis = read('/Users/william/Documents/lb01/14_328z.mseed')
            M=[[2,50,57],[8,22,15],[13,11,9],[19,7,58],[22,52,31]]
        if x==3:
            seis = read('/Users/william/Documents/lb01/14_329z.mseed')
            M=[[1,49,8],[2,25,56],[7,38,37],[8,56,5],[13,40,50]]
        if x==4:
            seis = read('/Users/william/Documents/lb01/14_330z.mseed')
            M=[[2,22,9],[3,14,20],[7,6,43],[11,6,13],[17,36,1]]
        if x==5:
            seis = read('/Users/william/Documents/lb01/14_331z.mseed')
            M=[[1,25,18],[12,59,20],[15,5,46],[19,11,43],[21,33,7]]
        if x==6:
            seis = read('/Users/william/Documents/lb01/14_332z.mseed')
            M=[[3,21,27],[4,20,18],[7,7,38],[12,42,54],[22,16,3]]
        if x==7:
            seis = read('/Users/william/Documents/lb01/14_333z.mseed')
            M=[[3,23,12],[4,6,44],[5,49,59],[15,51,20],[19,46,30]]
        if x==8:
            seis = read('/Users/william/Documents/lb01/14_334z.mseed') 
            M=[[1,9,2],[5,30,42],[12,31,5],[20,15,40],[22,27,20]]
        if x==9:
            seis = read('/Users/william/Documents/lb01/14_335z.mseed')
            M=[[1,22,40],[5,27,21],[10,13,51],[12,41,20],[22,1,41]]
        for i in range(0,len(M),1):   
            
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]
            window=65
            #
            seis[0].filter("bandpass", freqmin=fin,freqmax=fmax)
            start=seis[0].stats.starttime + h*60*60 + m*60 + s -10 # define start time (hours mins seconds) 
            end = start + window # define length of window in seconds
            trc = seis[0].slice(starttime = start  , endtime= end)  
            stream.append(trc)
            
    #%%        
    
    for x in range(0,50):
        trc=stream[x].normalize()
        stream2.append(trc)
    stack_norm = np.sum([abs(trc.data) for trc in stream2], axis=0)
    stack_norm = stack_norm/50
    
    return(stack_norm)