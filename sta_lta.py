#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 15:02:54 2018

@author: william
"""

from obspy.core import read
from obspy.signal.cross_correlation import correlate
import glob
import time
import obspy
import numpy as np
from numpy import argmax
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset

   #%% constants    
shift=100
step=2
#s_window=7*60*60 
window=40

#%% Read in waveforms      
sample = read('/Users/william/Documents/lb01/14_339z.mseed')
trace=sample[0]
trace.detrend(type='demean')
trace.filter("bandpass", freqmin=3,freqmax=10)

#window endpoints
start= sample[0].stats.starttime + 21*60*60 + 10*60  #+18*60*60 +20*60  #time window start 
end= start + 70*60
#end=sample[0].stats.endtime 
trs = trace.slice(starttime = start  , endtime= end) #cut out sample waveform with same window length as chosen event
#trs_e = obspy.signal.filter.envelope(trs.data)
#print('reference waveform')


trs.plot(type='relative',color='b')#, starttime=start , endtime=end)

sr = trace.stats.sampling_rate
nsta=int(2*sr)                                      #2
nlta=int(10*sr)                                     #20
stream=trs.data
cft=recursive_sta_lta(stream, nsta, nlta)
trig_on=6                                            #8
trig_off=0.2                                        #0.2
plot_trigger(trs, cft, trig_on, trig_off) 

on_off = trigger_onset(cft,trig_on,trig_off)

for x in range(0,len(on_off)):
    tr = trace.slice(starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
    tr_e = obspy.signal.filter.envelope(tr.data)
    

    #%% frequency info
    tr_data=tr.data
    m=np.mean(tr_data)
    tr_data = tr_data-m
    famp = abs(np.fft.fft(tr_data))
    
    # peak f
    peak= argmax(abs(famp))/window
    if peak > 50:
        peak = 100-peak

    # median f
    hal =sum(famp)/4 #half of the first half of the complete fft (ie. less than nq)
    num=0
    misf_i=1000000000
    for t in range(0,len(famp),1):
        n=famp[t]                  
        num += n
        misf = abs(num-hal)
        if misf < misf_i:
            misf_i = misf
            cf=t/window
    
    
#%% plot if not EXP  
    crite=0.25
#    if abs(top_v) < crite: #only allow waveforms that are not EXPs
    if 5 < cf :    
#        print(cf)
        trs.plot(type='relative',color='b', starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)









# DO NOT DELETE
