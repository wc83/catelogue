#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 14:26:07 2018

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
shift=200
step=2
crite_eq=0.10
crite=0.60
crite_n=0.10
#s_window=7*60*60 
window=30
fmin=3
fmax=12
file_sta = open('/Users/william/Documents/scanner/full_eq_scan_v1b.txt','w')
#%% Read in waveforms   
## Noise1
sample = read('/Users/william/Documents/lb01/14_336z.mseed')
trace=sample[0]
trace.detrend(type='demean')
trace.filter("bandpass", freqmin=fmin,freqmax=fmax)

startw= sample[0].stats.starttime + 14*60*60 + 11*60 + 34 #time window start - several s before Pwave arrival
endw=startw + window
trc = trace.slice(starttime = startw  , endtime= endw) #cut out sample waveform with same window length as chosen event
trc_e = obspy.signal.filter.envelope(trc.data)
#print('reference Noise(1) waveform')
#trc.plot(type='relative',color='b', starttime=startw , endtime=endw)

#EQ
trs2_e=stacked_eq(fmin,fmax)

#noise2
sample3 = read('/Users/william/Documents/lb01/14_335z.mseed')
trace3=sample3[0]
trace3.detrend(type='demean')
trace3.filter("bandpass", freqmin=fmin,freqmax=fmax)

start3= sample3[0].stats.starttime + 20*60*60 + 57*60 + 50 #time window start - several s before Pwave arrival
end3=start3 + window
trs3 = trace3.slice(starttime = start3  , endtime= end3) #cut out sample waveform with same window length as chosen event
trs3_e = obspy.signal.filter.envelope(trs3.data)
#print('reference noise(2) waveform')
#trs3.plot(type='relative',color='b', starttime=start3 , endtime=end3)
#%%  Read in days to scan
stream = sample.copy()
stream = stream.clear()
stream3 = sample.copy()
stream3 = stream3.clear()
for x in range(0,20):
    seis1=get_LB01z(x)
    seis3=get_LB03z(x)
    tr =seis1[0]
    tr3=seis3[0]
    tr.filter("bandpass", freqmin=fmin,freqmax=fmax)
    tr3.filter("bandpass", freqmin=fmin,freqmax=fmax)
    stream.append(tr)
    stream3.append(tr3)
#    print(x)
#stream.sort(keys=['starttime']) 
print('files read in') 
#print(stream)
#print(stream2)

#%%
sr = seis1[0].stats.sampling_rate
nsta=int(1*sr)                                      #2
nlta=int(20*sr)                                     #20
trig_on=7                                            #8
trig_off=0.2                                         #0.2

event=[]
event_other=[]
event_eq=[]
count=0
for t in range(0,len(stream),1): 
    print('Day',t+1,'of',len(stream))
    file_sta.write("\n")
    trace1=stream[t]
    trace3=stream3[t]

    
    #window endpoints
    start= stream[t].stats.starttime   #time window start 
 
#    trs.plot(type='relative',color='b')#, starttime=start , endtime=end)
    st=trace1.data
    st3=trace3.data
    cft=recursive_sta_lta(st, nsta, nlta)
    cft3=recursive_sta_lta(st3, nsta, nlta)
                                       
#    plot_trigger(trs2, cft2, trig_on, trig_off) 
    
    on_off1 = trigger_onset(cft,trig_on,trig_off)
    on_off3= trigger_onset(cft3,trig_on,trig_off)








#%% All LB01 triggered events
    event_day=[]
    for x in range(0,len(on_off1)):
        tr = trace1.slice(starttime=start+(on_off1[x,0]/sr) , endtime=start+(on_off1[x,1]/sr))
        tr_len=tr.stats.endtime - tr.stats.starttime 
    
        #%% frequency info

        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
        
        tr.detrend(type='demean')
        amp=abs(tr.max())
        whole = sum(abs(tr.data))
        av_amp = whole/(100*tr_len)
        amp_rat = amp/av_amp
        # only look for things above noise level 
        
        if amp > 300:
            if amp_rat < 18:
                if 4 < cf < 10 :  
                    if 3 < peak < 15:
                        if 3 < bwid < 14:
                            
#                            tr.plot(type='relative',color='g')
                            trs_e = obspy.signal.filter.envelope(tr.data)
                            top_v,top,corell = corel(trc_e,trs_e,shift)
                          
                            top_v_eq,top_eq,corell_eq = corel(trs2_e,trs_e,shift)
                            
                            top_v_n,top_n,corell_n = corel(trs3_e,trs_e,shift)
                       
                            if abs(top_v) < crite:
                                if abs(top_v_n) > crite_n:
                                    if abs(top_v_eq) > crite_eq:
#                                            print('EQ:')
#                                        tr.plot(type='relative',color='g')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                                        event_day.append(tr.stats.starttime)


#                                        file_sta.write(str(tr.stats.starttime + 8))
#                                        file_sta.write("\n")
#                                        event_eq.append(tr.stats.starttime) 
#                                        count+=1

    
#print(event)    
#file_sta.close()
#print(count)
               

    #%% LB03 only events
     
    LB03_start=[]
    LB03_end=[]                         
    for q in range(0,len(on_off3)):
        value = on_off3[q,0]
        
        near,ind=find_nearest(on_off1[:,0], value)
#        print('near =',near, abs(near-value))
        if abs(value-near) > 100*15: #15s          
#            print('seen in both',value)
            LB03_start.append(on_off3[q,0]/sr)
            LB03_end.append(on_off3[q,1]/sr)
#            print('lb03 finds =',on_off2[q,0]/sr)
    for x in range(0,len(LB03_start)):
        tr = trace3.slice(starttime=start+(LB03_start[x]) , endtime=start+(LB03_end[x]))
        tr.detrend(type='demean')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        
        #%% frequency info
        try:
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      

        
            amp=abs(tr.max())
            whole = sum(abs(tr.data))
            av_amp = whole/(100*tr_len)
            amp_rat = amp/av_amp
            # only look for things above noise level 
            
            if amp > 300:
                if amp_rat < 18:
                    if 4 < cf < 10 :  
                        if 3 < peak < 15:
                            if 3 < bwid < 14:
    #                            tr.plot(type='relative',color='g')
                                trs_e = obspy.signal.filter.envelope(tr.data)
                                top_v,top,corell = corel(trc_e,trs_e,shift)
                              
                                top_v_eq,top_eq,corell_eq = corel(trs2_e,trs_e,shift)
                                
                                top_v_n,top_n,corell_n = corel(trs3_e,trs_e,shift)
                           
                                if abs(top_v) < crite:
                                    if abs(top_v_n) > crite_n:
                                        if abs(top_v_eq) > crite_eq:
    #                                            print('EQ:')
    #                                        print('cf',cf,'peak',peak,'bwid',bwid,'Xc',top_v_eq,'rat',amp_rat)
    #                                        tr.plot(type='relative',color='g')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                                            event_day.append(tr.stats.starttime)
#                                            print('lb03 event added')
        except:
            print('find at edge of data')
    event_day.sort()
#    print(event_day)
    for p in range(0,len(event_day)):
        file_sta.write(str(event_day[p]))
        file_sta.write("\n")
file_sta.close()
print('end of scan')
        

        

    
