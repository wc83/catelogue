#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 16:22:49 2018

@author: william
"""

from obspy.core import read
from obspy.signal.cross_correlation import correlate
import time
import obspy
import numpy as np
from numpy import argmax
import matplotlib.pyplot as plt

T1= time.clock()
file_exp = open('/Users/william/Documents/scanner/full_exp_scan_v1c.txt','w')
   #%% constants    
shift=100
step=2
crite=0.55 
window=45
fmin=0.1
fmax=10
#%% Reference waveforms     
sample = read('/Users/william/Documents/lb01/14_330z.mseed')
trace=sample[0]
trace.detrend(type='demean')
trace.filter("bandpass", freqmin=fmin,freqmax=fmax)
start= sample[0].stats.starttime + 2*60*60 + 22*60 + 12 #time window start - several s before Pwave arrival
end=start + 45
trs = trace.slice(starttime = start  , endtime= end) #cut out sample waveform with same window length as chosen event
trs_e = obspy.signal.filter.envelope(trs.data)


#************ GET WAVEFORM FOR LB03 ************

#%% open empty streams for each station
stream1 = sample.copy()
stream1 = stream1.clear()
stream2 = sample.copy()
stream2 = stream2.clear()


#%% Read in Each station
for x in range(0,10):
    seis=get_LB01z(x)
    seis2=get_LB03z(x)
    tr =seis[0]
    tr2=seis2[0]
    stream1.append(tr)
    stream2.append(tr2)
print('files read in') 
#print(stream1)
#print(stream2)

#%% Open matracies

Events_per_day=[]
Event_store=[]
All_events=[]



#%% Loop over whole time
for p in range(0,len(stream1)):
    file_exp.write("\n")
    print('Day',p+1,'of',len(stream1))
    Events_today=np.zeros(shape=(100,2))
    todays_event=[]
    Number_of_events_today=0
    find_count=0
    
    trace=stream1[p]
    trace2=stream2[p]

#%% Station 1 scan    
    if trace.stats.npts > 202:
        begin = stream1[p].stats.starttime
        add=0
        lastcorr=0         
        dl=stream1[p].stats.endtime - stream1[p].stats.starttime
            
        for t in range(0,trace.stats.npts, step):
            view_start = begin + t + add
            view_end = view_start + window
            if view_start > stream1[p].stats.endtime - window :
                break
    
            trc = stream1[p].slice(starttime = view_start  , endtime= view_end)
            trc.detrend(type='demean')
            trc.filter("bandpass", freqmin=fmin,freqmax=fmax)
            # only look for things above noise level 
            amp=abs(trc.max())
            
            if amp > 1600:                
                peak,cf,bwid,bwid25 = freq_info(trc.data,trc.stats.starttime,trc.stats.endtime)         
                if 0.75 < cf < 2.75: 
                    if 0.2 < peak < 2.5:        
                        if 0.2 < bwid < 3: 
                                 #correlate between data and earthquake envelopes 
                            trc_e = obspy.signal.filter.envelope(trc.data) 
                            top_v,top,corell = corel(trs_e,trc_e,shift)
                            
                            if 0.0 in corell : # in the case of missing data - which will cause function to crash
                                break            
 
                #            if the correlation is positive, exp found
                            if abs(top_v) > crite:
                                find_count += 1
                                expt = view_start - (top/100) +2
                                rt=expt.timestamp	
                                jd=expt.julday
                                yr=expt.year
                                mo=expt.month
                                da=expt.day
                                hr=expt.hour
                                mi=expt.minute
                                se=expt.second
                                ms=expt.microsecond
                                row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,x])
                                print('lb01',rt)
                                
                                todays_event.append(expt)
                                Event_store.append(row)
                                Events_today[find_count][0]=rt
                                Events_today[find_count][1]=x
                                All_events.append(expt)
                                Number_of_events_today += 1
                                
                                add += 60-step               #skip 60s to avoid double catch

            else:
                add += 38
                lastcorr=0     
    
 #%% Station 2 scan  
    if trace2.stats.npts > 202:
        print('lb03 scan')
        begin = stream2[p].stats.starttime
        add=0
        lastcorr=0         
        dl=stream2[p].stats.endtime - stream2[p].stats.starttime
            
        for t in range(0,trace2.stats.npts, step):
            view_start = begin + t + add
            view_end = view_start + window
            if view_start > stream2[p].stats.endtime - window :
                break
    
            trc2 = stream2[p].slice(starttime = view_start  , endtime= view_end)
            trc2.detrend(type='demean')
            trc2.filter("bandpass", freqmin=fmin,freqmax=fmax)
            # only look for things above noise level 
            amp=abs(trc.max())
            
            if amp > 1600:       
                peak2,cf2,bwid2,bwid252 = freq_info(trc2.data,trc2.stats.starttime,trc2.stats.endtime)      
                if 0.75 < cf2 < 2.75: 
                    if 0.2 < peak2 < 2.5:        
                        if 0.2 < bwid2 < 3:
                                 #correlate between data and earthquake envelopes 
                            trc2_e = obspy.signal.filter.envelope(trc2.data) 
                            top_v,top,corell = corel(trs_e,trc2_e,shift)
                            
                            if 0.0 in corell : # in the case of missing data - which will cause function to crash
                                break            
 
                #            if the correlation is positive, exp found
                            if abs(top_v) > crite:
                                print('LB03 find')
                                expt = view_start - (top/100) +2
                                rt=expt.timestamp	
                                near,ind=find_nearest(Events_today[:][0], rt)
                                if abs(rt-near) > 15: #15s   *** check if time stamp works in secods ****
                                    print('LB03 unique',rt)
                                    find_count += 1
                                    jd=expt.julday
                                    yr=expt.year
                                    mo=expt.month
                                    da=expt.day
                                    hr=expt.hour
                                    mi=expt.minute
                                    se=expt.second
                                    ms=expt.microsecond
                                    row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,x])
                                    
                                    
                                    todays_event.append(expt)
                                    Event_store.append(row)
                                    Events_today[find_count][0]=rt
                                    Events_today[find_count][1]=x
                                    All_events.append(expt)
                                    Number_of_events_today += 1
                                    
                                add += 60-step               #skip 60s to avoid double catch

            else:
                add += 38
                lastcorr=0 
        Events_per_day.append(Number_of_events_today)
        Events_today.sort()
        for l in range(0,len(todays_event)):
            file_exp.write(str(todays_event[l]))
            file_exp.write("\n")
        
#%%
print("End of Scan")
T2=time.clock()
elapsed= T2-T1
print('Time taken:', elapsed)
file_exp.close()            
    
   
