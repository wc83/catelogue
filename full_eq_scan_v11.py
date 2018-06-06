#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 16:31:26 2018

@author: william
"""

from obspy.core import read
from obspy.signal.cross_correlation import correlate
import glob
import time
import obspy
import numpy as np
from numpy import argmax
from obspy import Stream
from obspy import UTCDateTime
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
T1=time.clock()

#%% read in explosions to remove
ts=np.zeros(shape=(1,1))
counter=0
with open('/Users/william/Documents/scanner/all_stations/EXP_coincidence_V1.csv') as f:
    reader = csv.reader(f)
    next(reader, None)
    for row in reader:
        ts[counter][0]=float(row[0])
        ts=np.lib.pad(ts, ((0,1),(0,0)), 'constant', constant_values=(0))
        counter+= 1


   #%% constants  
shift=500               # was 200 but increased to 500 to get higher Xcorr - will be slower
crite=0.5               # V8 =0.55 V8b = 0.6   
fmin=2
fmax=18
#%%
#file_sta = open('/Users/william/Documents/scanner/all_stations/full_eq_scan_v10e.txt','w')
file_cor2 = open('/Users/william/Documents/scanner/all_stations/full_correlated2_eq_scan_v11b.txt','w')
file_lb01= open('/Users/william/Documents/scanner/LB01/LB01_full_eq_scan_v11b.txt','w')
file_lb02= open('/Users/william/Documents/scanner/LB02/LB02_full_eq_scan_v11b.txt','w')
file_lb03= open('/Users/william/Documents/scanner/LB03/LB03_full_eq_scan_v11b.txt','w')
file_lb04= open('/Users/william/Documents/scanner/LB04/LB04_full_eq_scan_v11b.txt','w')
file_lb05= open('/Users/william/Documents/scanner/LB05/LB05_full_eq_scan_v11b.txt','w')
file_lb06= open('/Users/william/Documents/scanner/LB06/LB06_full_eq_scan_v11b.txt','w')

#%% Read in waveforms 
  
trc2_e1,trc2_e2,trc2_e3,trc2_e4,trc2_e5,trc2_e6=stacked_LBz_eq(fmin,fmax) 

#%% open empty streams for each station

stream1 = Stream()
stream2 = Stream()
stream3 = Stream()
stream4 = Stream()
stream5 = Stream()
stream6 = Stream()

#%% Read in Each station

num_active=[]

for x in range(14,17):
    
    seis1,seis2,seis3,seis4,seis5,seis6,num=get_LBz(x)    
    num_active.append(num)

    tr1=seis1[0]
    tr2=seis2[0]
    tr3=seis3[0]
    tr4=seis4[0]
    tr5=seis5[0]
    tr6=seis6[0]

    stream1.append(tr1)
    stream2.append(tr2)
    stream3.append(tr3)
    stream4.append(tr4)
    stream5.append(tr5)
    stream6.append(tr6)
   
print('files read in') 

#%%
sr = seis1[0].stats.sampling_rate
nsta=int(1*sr)                                      #1
nlta=int(20*sr)                                     #20
trig_on=11                                            #8 // 12 in 10d
trig_off=0.2                                         #0.2

count=0
corr_event2=np.zeros(shape=(1,15))
corr_event_count2=0


#%%

for t in range(0,len(stream1),1): 
    event_count=0

    print('Day',t+1,'of',len(stream1))
#    file_sta.write("\n")
    file_cor2.write("\n")
    corr_day2=[]
    
    trace1=stream1[t]
    trace2=stream2[t]
    trace3=stream3[t]
    trace4=stream4[t]
    trace5=stream5[t]
    trace6=stream6[t]
 
    st1=trace1.data
    st2=trace2.data
    st3=trace3.data
    st4=trace4.data
    st5=trace5.data
    st6=trace6.data
    
    cft1=recursive_sta_lta(st1, nsta, nlta)
    cft2=recursive_sta_lta(st2, nsta, nlta)
    cft3=recursive_sta_lta(st3, nsta, nlta)
    cft4=recursive_sta_lta(st4, nsta, nlta)
    cft5=recursive_sta_lta(st5, nsta, nlta)
    cft6=recursive_sta_lta(st6, nsta, nlta)
                                   
#    plot_trigger(trs2, cft2, trig_on, trig_off) 
    on_off_all=np.zeros(shape=(1,8))
    
    on_off1 = trigger_onset(cft1,trig_on,trig_off)
    on_off2 = trigger_onset(cft2,trig_on,trig_off)
    on_off3 = trigger_onset(cft3,trig_on,trig_off)
    on_off4 = trigger_onset(cft4,trig_on,trig_off)
    on_off5 = trigger_onset(cft5,trig_on,trig_off)
    on_off6 = trigger_onset(cft6,trig_on,trig_off)
    
#%% All LB01 triggered events
    event_day=[]

    for x in range(0,len(on_off1)):
        start=trace1.stats.starttime
        trace1.detrend(type='demean')
        trace1.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset1=start+(on_off1[x,0]/sr) 
        offset1= start+(on_off1[x,1]/sr) 
        event_len = offset1-onset1
        tr = trace1.slice(starttime=onset1 - 15 , endtime=offset1 +25)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
        
        
        if (4 < cf < 11) and (2 < peak < 12) and (2 < bwid < 12) and (4 < bwid25) and(event_len > 18) :
            rt=onset1.timestamp
            near,ind=find_nearest(ts[:,0], rt)
            if abs(near-rt) > 20:
                trs_e = obspy.signal.filter.envelope(tr.data)
                top_v_eq,top_eq,corell_eq = corel(trc2_e1,trs_e,shift) 
                 
                if abs(top_v_eq) > crite:
                
                    file_lb01.write(str(onset1))
                    file_lb01.write("\n")                                            
                    event_day.append(onset1)
                	
                    jd=onset1.julday
                    yr=onset1.year
                    mo=onset1.month
                    da=onset1.day
                    hr=onset1.hour
                    mi=onset1.minute
                    se=onset1.second
                    ms=onset1.microsecond/10000
                    
                    on_off_all[event_count][0]=rt
                    on_off_all[event_count][1]=event_len
                    on_off_all[event_count][2]=1
                    event_count += 1 
                    on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0)) 
                    
                    if num_active[t]==1:
                        if rt >0:
                            corr_event2[corr_event_count2][0]=rt
                            corr_event2[corr_event_count2][1]=event_len
                            corr_event2[corr_event_count2][2]=t
                            corr_event2[corr_event_count2][3]=yr
                            corr_event2[corr_event_count2][4]=jd
                            corr_event2[corr_event_count2][5]=hr
                            corr_event2[corr_event_count2][6]=mi
                            corr_event2[corr_event_count2][7]=se
                            corr_event2[corr_event_count2][8]=ms
                            corr_event2[corr_event_count2][9]=1                          
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset1)
                            
    #%% LB02 only events

    for x in range(0,len(on_off2)):
        start2=trace2.stats.starttime
        trace2.detrend(type='demean')
        trace2.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset2=start2+(on_off2[x,0]/sr) 
        offset2= start2+(on_off2[x,1]/sr)
        event_len = offset2-onset2
        tr = trace2.slice(starttime=onset2 -15 , endtime=offset2 +25)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      

        if (4 < cf < 11) and (2 < peak < 13) and (1.2 < bwid < 11) and (4 < bwid25 )  and (event_len > 18) :# and (ev_amp2 < amp) and (event_len > 18) :
            rt=onset2.timestamp	
            near,ind=find_nearest(ts[:,0], rt)
            if abs(near-rt) > 20: 
                   
                trs_e = obspy.signal.filter.envelope(tr.data)
                top_v_eq,top_eq,corell_eq = corel(trc2_e2,trs_e,shift)
    
                if abs(top_v_eq) > crite:
    
                    file_lb02.write(str(onset2))
                    file_lb02.write("\n")
                    
                    jd=onset2.julday
                    yr=onset2.year
                    mo=onset2.month
                    da=onset2.day
                    hr=onset2.hour
                    mi=onset2.minute
                    se=onset2.second
                    ms=onset2.microsecond/10000                                            
    
                    value = rt 
                    near,ind=find_nearest(on_off_all[:,0], value)
            
                    if abs(value-near) > 30: #if outside 30s of an event, save to full list  
                        on_off_all[event_count][0]=rt
                        on_off_all[event_count][1]=event_len
                        on_off_all[event_count][3]=1
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                        event_day.append(onset2)
                    else: 
                        on_off_all[ind][3]=1
                        
                    if num_active[t] > 1:    
                        if abs(value-near) < 10: # if already found, save as a 2x coincidence event
                            near_co,indx=find_nearest(corr_event2[:,0], value)
                            if abs(value-near_co)>30:  
                                if rt >0:  
                                    corr_event2[corr_event_count2][0]=rt
                                    corr_event2[corr_event_count2][1]=event_len
                                    corr_event2[corr_event_count2][2]=t
                                    corr_event2[corr_event_count2][3]=yr
                                    corr_event2[corr_event_count2][4]=jd
                                    corr_event2[corr_event_count2][5]=hr
                                    corr_event2[corr_event_count2][6]=mi
                                    corr_event2[corr_event_count2][7]=se
                                    corr_event2[corr_event_count2][8]=ms
                                    corr_event2[corr_event_count2][9]=on_off_all[ind][2]
                                    corr_event2[corr_event_count2][10]=on_off_all[ind][3] 
                                    corr_event_count2 += 1
                                    corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    corr_day2.append(onset2)
                    
                    else: 
                        if rt >0:
                            corr_event2[corr_event_count2][0]=rt
                            corr_event2[corr_event_count2][1]=event_len
                            corr_event2[corr_event_count2][2]=t
                            corr_event2[corr_event_count2][3]=yr
                            corr_event2[corr_event_count2][4]=jd
                            corr_event2[corr_event_count2][5]=hr
                            corr_event2[corr_event_count2][6]=mi
                            corr_event2[corr_event_count2][7]=se
                            corr_event2[corr_event_count2][8]=ms
                            corr_event2[corr_event_count2][10]=1   
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset2) 

#%% LB03 only events

    for x in range(0,len(on_off3)):
        start3=trace3.stats.starttime
        trace3.detrend(type='demean')
        trace3.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset3=start3+(on_off3[x,0]/sr) 
        offset3= start3+(on_off3[x,1]/sr)
        event_len = offset3-onset3
        tr = trace3.slice(starttime=onset3 -15, endtime=offset3 +25 )
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
    
        if (4 < cf < 13) and (3 < peak < 15) and (1.2 < bwid < 11) and (4 < bwid25 )  and (event_len > 18) :# and (ev_amp3 < amp) and (event_len > 18) :
            rt=onset3.timestamp	
            near,ind=find_nearest(ts[:,0], rt)
            if abs(near-rt) > 20:
                trs_e = obspy.signal.filter.envelope(tr.data)          
                top_v_eq,top_eq,corell_eq = corel(trc2_e3,trs_e,shift)            
           
                if abs(top_v_eq) > crite:
                    
                    file_lb03.write(str(onset3))
                    file_lb03.write("\n") 
                    
                    jd=onset3.julday
                    yr=onset3.year
                    mo=onset3.month
                    da=onset3.day
                    hr=onset3.hour
                    mi=onset3.minute
                    se=onset3.second
                    ms=onset3.microsecond/10000
                    
                    value = rt
                    near,ind=find_nearest(on_off_all[:,0], value)
    
                    if abs(value-near) > 30: #if outside 30s of an event, save to full list  
                        on_off_all[event_count][0]=rt
                        on_off_all[event_count][1]=event_len
                        on_off_all[event_count][4]=1
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                        event_day.append(onset3)
                    else: 
                        on_off_all[ind][4]=1
                        
                    if num_active[t] > 1:    
                        if abs(value-near) < 10: # if already found, save as a 2x coincidence event
                            near_co,indx=find_nearest(corr_event2[:,0], value)
                            if abs(value-near_co)>30:
                                if rt >0:    
                                    corr_event2[corr_event_count2][0]=rt
                                    corr_event2[corr_event_count2][1]=event_len
                                    corr_event2[corr_event_count2][2]=t
                                    corr_event2[corr_event_count2][3]=yr
                                    corr_event2[corr_event_count2][4]=jd
                                    corr_event2[corr_event_count2][5]=hr
                                    corr_event2[corr_event_count2][6]=mi
                                    corr_event2[corr_event_count2][7]=se
                                    corr_event2[corr_event_count2][8]=ms
                                    corr_event2[corr_event_count2][9]=on_off_all[ind][2]
                                    corr_event2[corr_event_count2][10]=on_off_all[ind][3]
                                    corr_event2[corr_event_count2][11]=on_off_all[ind][4]
                                    corr_event_count2 += 1
                                    corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    corr_day2.append(onset3) 
                            else:
                                corr_event2[indx][11]=1                        
                    else: 
                        if rt >0:
                            corr_event2[corr_event_count2][0]=rt
                            corr_event2[corr_event_count2][1]=event_len
                            corr_event2[corr_event_count2][2]=t
                            corr_event2[corr_event_count2][3]=yr
                            corr_event2[corr_event_count2][4]=jd
                            corr_event2[corr_event_count2][5]=hr
                            corr_event2[corr_event_count2][6]=mi
                            corr_event2[corr_event_count2][7]=se
                            corr_event2[corr_event_count2][8]=ms
                            corr_event2[corr_event_count2][11]=1   
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset3) 
                            
    #%% LB04 only events

    for x in range(0,len(on_off4)):
        start4=trace4.stats.starttime
        trace4.detrend(type='demean')
        trace4.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset4=start4+(on_off4[x,0]/sr) 
        offset4= start4+(on_off4[x,1]/sr) 
        event_len = offset4-onset4
        tr = trace4.slice(starttime=onset4 -15, endtime=offset4 +25)        
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
    
        if (3 < cf < 11) and (2 < peak < 12) and (2 < bwid < 13) and (4 < bwid25 )  and (event_len > 18) :# and (ev_amp4 < amp) and (event_len > 18) :
            rt=onset4.timestamp	
            near,ind=find_nearest(ts[:,0], rt)
            if abs(near-rt) > 20:                
                trs_e = obspy.signal.filter.envelope(tr.data)
                top_v_eq,top_eq,corell_eq = corel(trc2_e4,trs_e,shift)
                               
                if abs(top_v_eq) > crite:
                            
                    file_lb04.write(str(onset4))
                    file_lb04.write("\n") 
                    	
                    jd=onset4.julday
                    yr=onset4.year
                    mo=onset4.month
                    da=onset4.day
                    hr=onset4.hour
                    mi=onset4.minute
                    se=onset4.second
                    ms=onset4.microsecond/10000
    
                    value = rt
                    near,ind=find_nearest(on_off_all[:,0], value)
    #                    corr_near,indx=find_nearest(corr_event2[:,0], value)
    
                    if abs(value-near) > 30: #if outside 30s of an event, save to full list  
                        on_off_all[event_count][0]=rt
                        on_off_all[event_count][1]=event_len
                        on_off_all[event_count][5]=1
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                        event_day.append(onset4)
                    else: 
                        on_off_all[ind][5]=1
                        
                    if num_active[t] > 1:    
                        if abs(value-near) < 10: # if already found, save as a 2x coincidence event
                            near_co,indx=find_nearest(corr_event2[:,0], value)
                            if abs(value-near_co)>30:      
                                if rt >0:
                                    corr_event2[corr_event_count2][0]=rt
                                    corr_event2[corr_event_count2][1]=event_len
                                    corr_event2[corr_event_count2][2]=t
                                    corr_event2[corr_event_count2][3]=yr
                                    corr_event2[corr_event_count2][4]=jd
                                    corr_event2[corr_event_count2][5]=hr
                                    corr_event2[corr_event_count2][6]=mi
                                    corr_event2[corr_event_count2][7]=se
                                    corr_event2[corr_event_count2][8]=ms
                                    corr_event2[corr_event_count2][9]=on_off_all[ind][2]
                                    corr_event2[corr_event_count2][10]=on_off_all[ind][3]
                                    corr_event2[corr_event_count2][11]=on_off_all[ind][4]
                                    corr_event2[corr_event_count2][12]=on_off_all[ind][5]
                                    corr_event_count2 += 1
                                    corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    corr_day2.append(onset4) 
                            else:
                                corr_event2[indx][12]=1                        
                    else: 
                        if rt > 0:
                            corr_event2[corr_event_count2][0]=rt
                            corr_event2[corr_event_count2][1]=event_len
                            corr_event2[corr_event_count2][2]=t
                            corr_event2[corr_event_count2][3]=yr
                            corr_event2[corr_event_count2][4]=jd
                            corr_event2[corr_event_count2][5]=hr
                            corr_event2[corr_event_count2][6]=mi
                            corr_event2[corr_event_count2][7]=se
                            corr_event2[corr_event_count2][8]=ms
                            corr_event2[corr_event_count2][12]=1   
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset4)

    #%% LB05 only events

    for x in range(0,len(on_off5)):
        start5=trace5.stats.starttime
        trace5.detrend(type='demean')
        trace5.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset5=start5+(on_off5[x,0]/sr) 
        offset5= start5+(on_off5[x,1]/sr) 
        event_len = offset5-onset5
        tr = trace5.slice(starttime=onset5 -15, endtime=offset5 +25)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
 
        if (4 < cf < 11) and (2.5 < peak < 11) and (1.5 < bwid < 12) and (2.5 < bwid25 )  and (event_len > 18) :# and (ev_amp5 < amp) and (event_len > 18) :
            rt=onset5.timestamp	
            near,ind=find_nearest(ts[:,0], rt)
            if abs(near-rt) > 20:
                trs_e = obspy.signal.filter.envelope(tr.data)
                top_v_eq,top_eq,corell_eq = corel(trc2_e5,trs_e,shift)
           
                if abs(top_v_eq) > crite:
                    
                    file_lb05.write(str(onset5))
                    file_lb05.write("\n") 
                    
                    jd=onset5.julday
                    yr=onset5.year
                    mo=onset5.month
                    da=onset5.day
                    hr=onset5.hour
                    mi=onset5.minute
                    se=onset5.second
                    ms=onset5.microsecond/10000
                    
                    value = rt
                    near,ind=find_nearest(on_off_all[:,0], value)
    #                    corr_near,indx=find_nearest(corr_event2[:,0], value)
    
                    if abs(value-near) > 30: #if outside 30s of an event, save to full list  
                        on_off_all[event_count][0]=rt
                        on_off_all[event_count][1]=event_len
                        on_off_all[event_count][6]=1
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                        event_day.append(onset5)
                    else: 
                        on_off_all[ind][6]=1
                        
                    if num_active[t] > 1:    
                        if abs(value-near) < 10: # if already found, save as a 2x coincidence event
                            near_co,indx=find_nearest(corr_event2[:,0], value)
                            if abs(value-near_co)>30:   
                                if rt >0:
                                    corr_event2[corr_event_count2][0]=rt
                                    corr_event2[corr_event_count2][1]=event_len
                                    corr_event2[corr_event_count2][2]=t
                                    corr_event2[corr_event_count2][3]=yr
                                    corr_event2[corr_event_count2][4]=jd
                                    corr_event2[corr_event_count2][5]=hr
                                    corr_event2[corr_event_count2][6]=mi
                                    corr_event2[corr_event_count2][7]=se
                                    corr_event2[corr_event_count2][8]=ms
                                    corr_event2[corr_event_count2][9]=on_off_all[ind][2]
                                    corr_event2[corr_event_count2][10]=on_off_all[ind][3]
                                    corr_event2[corr_event_count2][11]=on_off_all[ind][4]
                                    corr_event2[corr_event_count2][12]=on_off_all[ind][5]
                                    corr_event2[corr_event_count2][13]=on_off_all[ind][6]
                                    corr_event_count2 += 1
                                    corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    corr_day2.append(onset5)
                            else:
                                corr_event2[indx][13]=1
                    else:
                        if rt >0:
                            corr_event2[corr_event_count2][0]=rt
                            corr_event2[corr_event_count2][1]=event_len
                            corr_event2[corr_event_count2][2]=t
                            corr_event2[corr_event_count2][3]=yr
                            corr_event2[corr_event_count2][4]=jd
                            corr_event2[corr_event_count2][5]=hr
                            corr_event2[corr_event_count2][6]=mi
                            corr_event2[corr_event_count2][7]=se
                            corr_event2[corr_event_count2][8]=ms
                            corr_event2[corr_event_count2][13]=1   
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset5)
          
    #%% LB06 only events
                        
    for x in range(0,len(on_off6)):
        start6=trace6.stats.starttime
        trace6.detrend(type='demean')
        trace6.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset6=start6+(on_off6[x,0]/sr) 
        offset6= start6+(on_off6[x,1]/sr)
        event_len = offset6-onset6
        tr = trace6.slice(starttime=onset6 -15, endtime=offset6 +25)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
    
        if (4 < cf < 11) and (3 < peak < 11) and (2.5 < bwid < 9.5) and (4 < bwid25 )  and (event_len > 18) :# and (ev_amp5 < amp) and (event_len > 18) :  
            rt=onset6.timestamp	
            near,ind=find_nearest(ts[:,0], rt)
            
            if abs(near-rt) > 20:
                trs_e = obspy.signal.filter.envelope(tr.data)
                top_v_eq,top_eq,corell_eq = corel(trc2_e6,trs_e,shift)
                            
                if abs(top_v_eq) > crite:
                 
                    file_lb06.write(str(onset6))
                    file_lb06.write("\n") 
                    
                    jd=onset6.julday
                    yr=onset6.year
                    mo=onset6.month
                    da=onset6.day
                    hr=onset6.hour
                    mi=onset6.minute
                    se=onset6.second
                    ms=onset6.microsecond/10000
                    
                    value = rt
                    near,ind=find_nearest(on_off_all[:,0], value)
    #                    corr_near,indx=find_nearest(corr_event2[:,0], value)
    
                    if abs(value-near) > 30: #if outside 30s of an event, save to full list  
                        on_off_all[event_count][0]=rt
                        on_off_all[event_count][1]=event_len
                        on_off_all[event_count][7]=1
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                        event_day.append(onset6)
                    else: 
                        on_off_all[ind][7]=1
                        
                    if num_active[t] > 1:    
                        if abs(value-near) < 10: # if already found, save as a 2x coincidence event
                            near_co,indx=find_nearest(corr_event2[:,0], value)
                            if abs(value-near_co)>30:  
                                if rt > 0:                                  
                                    corr_event2[corr_event_count2][0]=rt
                                    corr_event2[corr_event_count2][1]=event_len
                                    corr_event2[corr_event_count2][2]=t
                                    corr_event2[corr_event_count2][3]=yr
                                    corr_event2[corr_event_count2][4]=jd
                                    corr_event2[corr_event_count2][5]=hr
                                    corr_event2[corr_event_count2][6]=mi
                                    corr_event2[corr_event_count2][7]=se
                                    corr_event2[corr_event_count2][8]=ms
                                    corr_event2[corr_event_count2][9]=on_off_all[ind][2]
                                    corr_event2[corr_event_count2][10]=on_off_all[ind][3]
                                    corr_event2[corr_event_count2][11]=on_off_all[ind][4]
                                    corr_event2[corr_event_count2][12]=on_off_all[ind][5]
                                    corr_event2[corr_event_count2][13]=on_off_all[ind][6]
                                    corr_event2[corr_event_count2][14]=on_off_all[ind][7]
                                    corr_event_count2 += 1
                                    corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    corr_day2.append(onset6)
                            else:
                                corr_event2[indx][14]=1
                    
                    else:
                        if rt > 0:
                            corr_event2[corr_event_count2][0]=rt
                            corr_event2[corr_event_count2][1]=event_len
                            corr_event2[corr_event_count2][2]=t
                            corr_event2[corr_event_count2][3]=yr
                            corr_event2[corr_event_count2][4]=jd
                            corr_event2[corr_event_count2][5]=hr
                            corr_event2[corr_event_count2][6]=mi
                            corr_event2[corr_event_count2][7]=se
                            corr_event2[corr_event_count2][8]=ms
                            corr_event2[corr_event_count2][14]=1   
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset6) 
       
#%%            
    event_day.sort()
    corr_day2.sort() 
    
    for p in range(0,len(corr_day2)):
        file_cor2.write(str(corr_day2[p]))
        file_cor2.write("\n")  
            

#%%             
col=0       
corr_event2=corr_event2[np.argsort(corr_event2[:,col])] 
if corr_event2[0][0]== 0.0:
    corr_event2 = np.delete(corr_event2, (0), axis=0)       
np.savetxt("/Users/william/Documents/scanner/all_stations/EQ_coincidence_V2b.csv", corr_event2,delimiter=",",header="Time_stamp,Event_length,Day_number,Year,Day_of_year,Hour,Min,Sec,Milisec,LBO1,LB02,LB03,LB04,LB05,LB06")

data_stream = Stream()
for r in range(0,len(corr_event2)):
    if corr_event2[r][9]==1:
        tr1a=stream1[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr1_save = tr1a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr1_save.detrend(type='demean')
        tr1_save.detrend(type='linear')
        data_stream.append(tr1_save)
    if corr_event2[r][10]==1:
        tr2a=stream2[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr2_save = tr2a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr2_save.detrend(type='demean')
        tr2_save.detrend(type='linear')
        data_stream.append(tr2_save)    
    if corr_event2[r][11]==1:
        tr3a=stream3[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr3_save = tr3a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr3_save.detrend(type='demean')
        tr3_save.detrend(type='linear')
        data_stream.append(tr3_save)         
    if corr_event2[r][12]==1:
        tr4a=stream4[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr4_save = tr4a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr4_save.detrend(type='demean')
        tr4_save.detrend(type='linear')
        data_stream.append(tr4_save)    
    if corr_event2[r][13]==1:
        tr5a=stream5[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr5_save = tr5a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr5_save.detrend(type='demean')
        tr5_save.detrend(type='linear')
        data_stream.append(tr5_save)
    if corr_event2[r][14]==1:
        tr6a=stream6[int(corr_event2[r][2])] 
        dstart=UTCDateTime(corr_event2[r][0])
        tr6_save = tr6a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr6_save.detrend(type='demean')
        tr6_save.detrend(type='linear')
        data_stream.append(tr6_save)

data_stream.write('/Users/william/Documents/scanner/output_data/EQ_data_stream_v2b.mseed', format='MSEED')         

        
#file_sta.close()
file_cor2 .close()
file_lb01.close()
file_lb02.close()
file_lb03.close()
file_lb04.close()
file_lb05.close()
file_lb06.close()

print('end of scan')
T2=time.clock()
elapsed= T2-T1
print('Time taken:', elapsed)       

        

    
