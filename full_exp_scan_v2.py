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
file_exp = open('/Users/william/Documents/scanner/full_exp_scan_v2.txt','w')
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
stream3 = sample.copy()
stream3 = stream3.clear()
stream4 = sample.copy()
stream4 = stream4.clear()
stream5 = sample.copy()
stream5 = stream5.clear()
stream6 = sample.copy()
stream6 = stream6.clear()
stream7 = sample.copy()
stream7 = stream7.clear()
stream8 = sample.copy()
stream8 = stream8.clear()

#%% Read in Each station
for x in range(0,25):
    seis1,seis2,seis3,seis4,seis5,seis6,seis7,seis8=get_LBz(x)

    
    tr1=seis1[0]
    tr2=seis2[0]
    tr3=seis3[0]
    tr4=seis4[0]
    tr5=seis5[0]
    tr6=seis6[0]
    tr7=seis7[0]
    tr8=seis8[0]
    
    stream1.append(tr1)
    stream2.append(tr2)
    stream3.append(tr3)
    stream4.append(tr4)
    stream5.append(tr5)
    stream6.append(tr6)
    stream7.append(tr7)
    stream8.append(tr8)
    
print('files read in') 

#%% Open matracies

Events_per_day=[]
Event_store=[]
All_events=[]


#%% Loop over whole time one day at a time
for p in range(0,len(stream1)):
    file_exp.write("\n")
    print('Day',p+1,'of',len(stream1))
    Events_today=np.zeros(shape=(1,2))
    todays_event=[]
    Number_of_events_today=0
    find_count=0
    
    trace1=stream1[p]
    trace2=stream2[p]
    trace3=stream3[p]
    trace4=stream4[p]
    trace5=stream5[p]
    trace6=stream6[p]
    trace7=stream7[p]

#%% Station LB01 scan    
    if trace1.stats.npts > 210:
        begin = trace1.stats.starttime
        add=0
        lastcorr=0         
        dl=trace1.stats.endtime - trace1.stats.starttime
            
        # Loop over single day, scanning for waveforms
        for t in range(0,trace1.stats.npts, step):
            view_start = begin + t + add
            view_end = view_start + window
            if view_start > trace1.stats.endtime - window :
                break
    
            trc = trace1.slice(starttime = view_start  , endtime= view_end)
            trc.detrend(type='linear')
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
                                
                                
                                todays_event.append(expt)
                                Event_store.append(row)
                                Events_today[find_count][0]=rt
                                Events_today[find_count][1]=x
                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                All_events.append(expt)
                                find_count += 1
                                Number_of_events_today += 1
                                
                                add += 60-step               #skip 60s to avoid double catch

            else:
                add += 38
                lastcorr=0     
    
 #%% Station 2 scan  
    if trace2.stats.npts > 210:
        begin = trace2.stats.starttime
        add=0
        lastcorr=0         
        dl=trace2.stats.endtime - trace2.stats.starttime
            
        for t in range(0,trace2.stats.npts, step):
            view_start = begin + t + add
            view_end = view_start + window
            if view_start > trace2.stats.endtime - window :
                break
    
            trc2 = trace2.slice(starttime = view_start  , endtime= view_end)
            trc2.detrend(type='linear')
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
                                expt = view_start - (top/100) +2
                                
                                rt=expt.timestamp	
                                near,ind=find_nearest(Events_today[:,0], rt)
                                if abs(rt-near) > 60: #60s
      
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
                                    
                                    todays_event.append(expt)
                                    Event_store.append(row)
                                    Events_today[find_count][0]=rt
                                    Events_today[find_count][1]=x
                                    Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    All_events.append(expt)
                                    find_count += 1
                                    Number_of_events_today += 1
                                    
                                add += 60-step               #skip 60s to avoid double catch

            else:
                add += 38
                lastcorr=0 
 
 #%% Station 3 scan  
    if trace3.stats.npts > 210:
        begin = trace3.stats.starttime
        add=0
        lastcorr=0         
        dl=trace3.stats.endtime - trace3.stats.starttime
            
        for t in range(0,trace3.stats.npts, step):
            view_start = begin + t + add
            view_end = view_start + window
            if view_start > trace3.stats.endtime - window :
                break
    
            trc3 = trace3.slice(starttime = view_start  , endtime= view_end)
            trc3.detrend(type='linear')
            trc3.detrend(type='demean')
            trc3.filter("bandpass", freqmin=fmin,freqmax=fmax)
            # only look for things above noise level 
            amp=abs(trc.max())
            
            if amp > 1600:       
                peak3,cf3,bwid3,bwid253 = freq_info(trc3.data,trc3.stats.starttime,trc3.stats.endtime)      
                if 0.75 < cf3 < 2.75: 
                    if 0.2 < peak3 < 2.5:        
                        if 0.2 < bwid3 < 3: 
                                 #correlate between data and earthquake envelopes 
                            trc3_e = obspy.signal.filter.envelope(trc3.data) 
                            top_v,top,corell = corel(trs_e,trc3_e,shift)
                            
                            if 0.0 in corell : # in the case of missing data - which will cause function to crash
                                break            
 
                #            if the correlation is positive, exp found
                            if abs(top_v) > crite:
                                expt = view_start - (top/100) +2
                                
                                rt=expt.timestamp	
                                near,ind=find_nearest(Events_today[:,0], rt)
                                if abs(rt-near) > 60: #60s
      
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
                                    
                                    todays_event.append(expt)
                                    Event_store.append(row)
                                    Events_today[find_count][0]=rt
                                    Events_today[find_count][1]=x
                                    Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    All_events.append(expt)
                                    find_count += 1
                                    Number_of_events_today += 1
                                    
                                add += 60-step               #skip 60s to avoid double catch

            else:
                add += 38
                lastcorr=0 
                
                 #%% Station 4 scan  
    if trace4.stats.npts > 210:
        begin = trace4.stats.starttime
        add=0
        lastcorr=0         
        dl=trace4.stats.endtime - trace4.stats.starttime
            
        for t in range(0,trace4.stats.npts, step):
            view_start = begin + t + add
            view_end = view_start + window
            if view_start > trace4.stats.endtime - window :
                break
    
            trc4 = trace4.slice(starttime = view_start  , endtime= view_end)
            trc4.detrend(type='linear')
            trc4.detrend(type='demean')
            trc4.filter("bandpass", freqmin=fmin,freqmax=fmax)
            # only look for things above noise level 
            amp=abs(trc.max())
            
            if amp > 1600:       
                peak4,cf4,bwid4,bwid254 = freq_info(trc4.data,trc4.stats.starttime,trc4.stats.endtime)      
                if 0.75 < cf4 < 2.75: 
                    if 0.2 < peak4 < 2.5:        
                        if 0.2 < bwid4 < 3: 
                                 #correlate between data and earthquake envelopes 
                            trc4_e = obspy.signal.filter.envelope(trc4.data) 
                            top_v,top,corell = corel(trs_e,trc4_e,shift)
                            
                            if 0.0 in corell : # in the case of missing data - which will cause function to crash
                                break            
 
                #            if the correlation is positive, exp found
                            if abs(top_v) > crite:
                                expt = view_start - (top/100) +2
                                
                                rt=expt.timestamp	
                                near,ind=find_nearest(Events_today[:,0], rt)
                                if abs(rt-near) > 60: #60s
      
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
                                    
                                    todays_event.append(expt)
                                    Event_store.append(row)
                                    Events_today[find_count][0]=rt
                                    Events_today[find_count][1]=x
                                    Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    All_events.append(expt)
                                    find_count += 1
                                    Number_of_events_today += 1
                                    
                                add += 60-step               #skip 60s to avoid double catch

            else:
                add += 38
                lastcorr=0 
                
    #%% Station 5 scan  
    if trace5.stats.npts > 210:
        begin = trace5.stats.starttime
        add=0
        lastcorr=0         
        dl=trace5.stats.endtime - trace5.stats.starttime
            
        for t in range(0,trace5.stats.npts, step):
            view_start = begin + t + add
            view_end = view_start + window
            if view_start > trace5.stats.endtime - window :
                break
    
            trc5 = trace5.slice(starttime = view_start  , endtime= view_end)
            trc5.detrend(type='linear')
            trc5.detrend(type='demean')
            trc5.filter("bandpass", freqmin=fmin,freqmax=fmax)
            # only look for things above noise level 
            amp=abs(trc.max())
            
            if amp > 1600:       
                peak5,cf5,bwid5,bwid255 = freq_info(trc5.data,trc5.stats.starttime,trc5.stats.endtime)      
                if 0.75 < cf5 < 2.75: 
                    if 0.2 < peak5 < 2.5:        
                        if 0.2 < bwid5 < 3: 
                                 #correlate between data and earthquake envelopes 
                            trc5_e = obspy.signal.filter.envelope(trc5.data) 
                            top_v,top,corell = corel(trs_e,trc5_e,shift)
                            
                            if 0.0 in corell : # in the case of missing data - which will cause function to crash
                                break            
 
                #            if the correlation is positive, exp found
                            if abs(top_v) > crite:
                                expt = view_start - (top/100) +2
                                
                                rt=expt.timestamp	
                                near,ind=find_nearest(Events_today[:,0], rt)
                                if abs(rt-near) > 60: #60s
      
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
                                    
                                    todays_event.append(expt)
                                    Event_store.append(row)
                                    Events_today[find_count][0]=rt
                                    Events_today[find_count][1]=x
                                    Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    All_events.append(expt)
                                    find_count += 1
                                    Number_of_events_today += 1
                                    
                                add += 60-step               #skip 60s to avoid double catch

            else:
                add += 38
                lastcorr=0 
                
    #%% Station 6 scan  
    if trace6.stats.npts > 210:
        begin = trace6.stats.starttime
        add=0
        lastcorr=0         
        dl=trace6.stats.endtime - trace6.stats.starttime
            
        for t in range(0,trace6.stats.npts, step):
            view_start = begin + t + add
            view_end = view_start + window
            if view_start > trace6.stats.endtime - window :
                break
    
            trc6 = trace6.slice(starttime = view_start  , endtime= view_end)
            trc6.detrend(type='linear')
            trc6.detrend(type='demean')
            trc6.filter("bandpass", freqmin=fmin,freqmax=fmax)
            # only look for things above noise level 
            amp=abs(trc.max())
            
            if amp > 1600:       
                peak6,cf6,bwid6,bwid256 = freq_info(trc6.data,trc6.stats.starttime,trc6.stats.endtime)      
                if 0.75 < cf6 < 2.75: 
                    if 0.2 < peak6 < 2.5:        
                        if 0.2 < bwid6 < 3: 
                                 #correlate between data and earthquake envelopes 
                            trc6_e = obspy.signal.filter.envelope(trc6.data) 
                            top_v,top,corell = corel(trs_e,trc6_e,shift)
                            
                            if 0.0 in corell : # in the case of missing data - which will cause function to crash
                                break            
 
                #            if the correlation is positive, exp found
                            if abs(top_v) > crite:
                                expt = view_start - (top/100) +2
                                
                                rt=expt.timestamp	
                                near,ind=find_nearest(Events_today[:,0], rt)
                                if abs(rt-near) > 60: #60s
      
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
                                    
                                    todays_event.append(expt)
                                    Event_store.append(row)
                                    Events_today[find_count][0]=rt
                                    Events_today[find_count][1]=x
                                    Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    All_events.append(expt)
                                    find_count += 1
                                    Number_of_events_today += 1
                                    
                                add += 60-step               #skip 60s to avoid double catch

            else:
                add += 38
                lastcorr=0 
                
    #%% Station 7 scan  
    if trace7.stats.npts > 210:
        begin = trace7.stats.starttime
        add=0
        lastcorr=0         
        dl=trace7.stats.endtime - trace7.stats.starttime
            
        for t in range(0,trace7.stats.npts, step):
            view_start = begin + t + add
            view_end = view_start + window
            if view_start > trace7.stats.endtime - window :
                break
    
            trc7 = trace7.slice(starttime = view_start  , endtime= view_end)
            trc7.detrend(type='linear')
            trc7.detrend(type='demean')
            trc7.filter("bandpass", freqmin=fmin,freqmax=fmax)
            # only look for things above noise level 
            amp=abs(trc.max())
            
            if amp > 1600:       
                peak7,cf7,bwid7,bwid257 = freq_info(trc7.data,trc7.stats.starttime,trc7.stats.endtime)      
                if 0.75 < cf7 < 2.75: 
                    if 0.2 < peak7 < 2.5:        
                        if 0.2 < bwid7 < 3: 
                                 #correlate between data and earthquake envelopes 
                            trc7_e = obspy.signal.filter.envelope(trc7.data) 
                            top_v,top,corell = corel(trs_e,trc7_e,shift)
                            
                            if 0.0 in corell : # in the case of missing data - which will cause function to crash
                                break            
 
                #            if the correlation is positive, exp found
                            if abs(top_v) > crite:
                                expt = view_start - (top/100) +2

                                rt=expt.timestamp	
                                near,ind=find_nearest(Events_today[:,0], rt)
                                if abs(rt-near) > 60: #60s
      
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
                                    
                                    todays_event.append(expt)
                                    Event_store.append(row)
                                    Events_today[find_count][0]=rt
                                    Events_today[find_count][1]=x
                                    Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    All_events.append(expt)
                                    find_count += 1
                                    Number_of_events_today += 1
                                    
                                add += 60-step               #skip 60s to avoid double catch

            else:
                add += 38
                lastcorr=0 
                
       #%% Station 8 scan  
    if trace8.stats.npts > 210:
        begin = trace8.stats.starttime
        add=0
        lastcorr=0         
        dl=trace8.stats.endtime - trace8.stats.starttime
            
        for t in range(0,trace8.stats.npts, step):
            view_start = begin + t + add
            view_end = view_start + window
            if view_start > trace7.stats.endtime - window :
                break
    
            trc8 = trace8.slice(starttime = view_start  , endtime= view_end)
            trc8.detrend(type='linear')
            trc8.detrend(type='demean')
            trc8.filter("bandpass", freqmin=fmin,freqmax=fmax)
            # only look for things above noise level 
            amp=abs(trc.max())
            
            if amp > 1600:       
                peak8,cf8,bwid8,bwid258 = freq_info(trc8.data,trc8.stats.starttime,trc8.stats.endtime)      
                if 0.75 < cf8 < 2.75: 
                    if 0.2 < peak8 < 2.5:        
                        if 0.2 < bwid8 < 3: 
                                 #correlate between data and earthquake envelopes 
                            trc8_e = obspy.signal.filter.envelope(trc8.data) 
                            top_v,top,corell = corel(trs_e,trc8_e,shift)
                            
                            if 0.0 in corell : # in the case of missing data - which will cause function to crash
                                break            
 
                #            if the correlation is positive, exp found
                            if abs(top_v) > crite:
                                expt = view_start - (top/100) +2

                                rt=expt.timestamp	
                                near,ind=find_nearest(Events_today[:,0], rt)
                                if abs(rt-near) > 60: #60s
      
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
                                    
                                    todays_event.append(expt)
                                    Event_store.append(row)
                                    Events_today[find_count][0]=rt
                                    Events_today[find_count][1]=x
                                    Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    All_events.append(expt)
                                    find_count += 1
                                    Number_of_events_today += 1
                                    
                                add += 60-step               #skip 60s to avoid double catch

            else:
                add += 38
                lastcorr=0 
                
                  
                
#%% collect information from each station for single day and save in arrays                
    Events_per_day.append(Number_of_events_today)
    todays_event.sort()
    for l in range(0,len(todays_event)):
        file_exp.write(str(todays_event[l]))
        file_exp.write("\n")
        
#%%
print("End of Scan")
T2=time.clock()
elapsed= T2-T1
print('Time taken:', elapsed)
file_exp.close()            
    
   
    
    
    
    
    
    
    
    