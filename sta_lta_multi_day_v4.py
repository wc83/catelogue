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
crit=0.4
crite=0.4
#s_window=7*60*60 
window=40
fmin=3
fmax=10
file_sta = open('/Users/william/Documents/scanner/sta_lta_EQ_V4.txt','w')
#%% Read in waveforms      
sample = read('/Users/william/Documents/lb01/14_330z.mseed')
stream = sample.copy()
stream = stream.clear()

for file in glob.glob('/Users/william/Documents/lb01/14_*.mseed'):

    seis = read(file)
    tr =seis[0]
    tr.detrend(type='demean')
    tr.filter("bandpass", freqmin=fmin,freqmax=fmax)
    stream.append(tr)
stream.sort(keys=['starttime']) 
print('files read in') 

sr = seis[0].stats.sampling_rate
nsta=int(1*sr)                                      #2
nlta=int(20*sr)                                     #20
trig_on=10                                            #8
trig_off=0.2                                         #0.2

event=[]
for t in range(0,len(stream),1): 
    print('Day',t+1,'of',len(stream))
    file_sta.write("\n")
    trs=stream[t]
    trace=trs
    
    #window endpoints
    start= stream[t].stats.starttime   #time window start 
#    end=stream[t].stats.endtime 
#    trs = trace.slice(starttime = start  , endtime= end) #cut out sample waveform with same window length as chosen event
    #trs_e = obspy.signal.filter.envelope(trs.data)
    #print('reference waveform')
    
    
#    trs.plot(type='relative',color='b')#, starttime=start , endtime=end)
    st=trs.data
    cft=recursive_sta_lta(st, nsta, nlta)
                                       
#    plot_trigger(trs, cft, trig_on, trig_off) 
    
    on_off = trigger_onset(cft,trig_on,trig_off)
    
    for x in range(0,len(on_off)):
        tr = trace.slice(starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
        
        tr_len=tr.stats.endtime - tr.stats.starttime 
    
        #%% frequency info
        tr_data=tr.data
        m=np.mean(tr_data)
        tr_data = tr_data-m
        famp = abs(np.fft.fft(tr_data))
        
        # peak f
        peak= argmax(abs(famp))/tr_len
        if peak > 50:
            peak = 100-peak
    
        # median f
        hal =sum(famp)/4 #half of the first half of the complete fft (ie. less than nq)
        num=0
        misf_i=1000000000
        for p in range(0,len(famp),1):
            n=famp[p]                  
            num += n
            misf = abs(num-hal)
            if misf < misf_i:
                misf_i = misf
                cf=p/tr_len
                
        
        tr.detrend(type='demean')
        amp=abs(tr.max())
        # only look for things above noise level 
        
        if amp > 300:
            if 5 < cf :  
#                print(cf, peak)
                if 12.4 > peak > 4.6:
                #plot 
#                    tr.plot(type='relative',color='b')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                    event.append(tr.stats.starttime)
                    file_sta.write(str(start+(on_off[x,0]/sr)))
                    file_sta.write("\n")
    
    
#print(event)    
file_sta.close()
    
    
    # DO NOT DELETE
