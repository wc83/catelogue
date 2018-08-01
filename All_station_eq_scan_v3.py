#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 11:27:07 2018

@author: william
"""


import time
import obspy
import numpy as np
from numpy import argmax
from obspy import Stream
from obspy import UTCDateTime
from obspy.signal.trigger import recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
import csv
T1=time.clock()

#%% read in explosions to remove
ts=np.zeros(shape=(1,1))
counter=0
with open('/Users/william/Documents/scanner/all_stations/EXP_all_coincidence_2_month_30.csv') as f:
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
file_cor2 = open('/Users/william/Documents/scanner/all_stations/All_correlated2_eq_scan_month_30.txt','w')
file_lb01= open('/Users/william/Documents/scanner/LB01/LB01_all_eq_scan_2_month_30.txt','w')
file_lb02= open('/Users/william/Documents/scanner/LB02/LB02_all_eq_scan_2_month_30.txt','w')
file_lb03= open('/Users/william/Documents/scanner/LB03/LB03_all_eq_scan_2_month_30.txt','w')
file_lb04= open('/Users/william/Documents/scanner/LB04/LB04_all_eq_scan_2_month_30.txt','w')
file_lb05= open('/Users/william/Documents/scanner/LB05/LB05_all_eq_scan_2_month_30.txt','w')
file_lb06= open('/Users/william/Documents/scanner/LB06/LB06_all_eq_scan_2_month_30.txt','w')
file_ls01= open('/Users/william/Documents/scanner/LS01/LS01_all_eq_scan_2_month_30.txt','w')
file_ls02= open('/Users/william/Documents/scanner/LS02/LS02_all_eq_scan_2_month_30.txt','w')
file_ls03= open('/Users/william/Documents/scanner/LS03/LS03_all_eq_scan_2_month_30.txt','w')
file_ls04= open('/Users/william/Documents/scanner/LS04/LS04_all_eq_scan_2_month_30.txt','w')
file_ls05= open('/Users/william/Documents/scanner/LS05/LS05_all_eq_scan_2_month_30.txt','w')
file_ls06= open('/Users/william/Documents/scanner/LS06/LS06_all_eq_scan_2_month_30.txt','w')

#%% Read in waveforms 
  
trc2_e1,trc2_e2,trc2_e3,trc2_e4,trc2_e5,trc2_e6=stacked_LBz_eq(fmin,fmax) 
trc2_s1,trc2_s2,trc2_s3,trc2_s4,trc2_s5,trc2_s6=stacked_LS_eq(fmin,fmax) 
print('Reference Waveforms Read In')

#%% open empty streams for each station

stream1 = Stream()
stream2 = Stream()
stream3 = Stream()
stream4 = Stream()
stream5 = Stream()
stream6 = Stream()
stream1s = Stream()
stream2s = Stream()
stream3s = Stream()
stream4s = Stream()
stream5s = Stream()
stream6s = Stream()

#%% Read in Each station

num_active=[]
lb_num_active=[]
for x in range(870,900):
    
    seis1,seis2,seis3,seis4,seis5,seis6,seis1s,seis2s,seis3s,seis4s,seis5s,seis6s,num,lb_num=get_all_stations(x)    
    num_active.append(num)
    lb_num_active.append(lb_num)
    
    tr1=seis1[0]
    tr2=seis2[0]
    tr3=seis3[0]
    tr4=seis4[0]
    tr5=seis5[0]
    tr6=seis6[0]
    tr1s=seis1s[0]
    tr2s=seis2s[0]
    tr3s=seis3s[0]
    tr4s=seis4s[0]
    tr5s=seis5s[0]
    tr6s=seis6s[0]

    stream1.append(tr1)
    stream2.append(tr2)
    stream3.append(tr3)
    stream4.append(tr4)
    stream5.append(tr5)
    stream6.append(tr6)
    stream1s.append(tr1s)
    stream2s.append(tr2s)
    stream3s.append(tr3s)
    stream4s.append(tr4s)
    stream5s.append(tr5s)
    stream6s.append(tr6s)
   
print('files read in') 
print('num_active =',num_active)
print('lb_num_active =',lb_num_active)
#%%
sr = seis1[0].stats.sampling_rate
nsta=int(1*sr)                                      #1
nlta=int(20*sr)                                     #20
trig_on=9                                          #10 in V2b, 9 in V2bb
trig_off=0.2                                         #0.2

count=0
corr_event2=np.zeros(shape=(1,21))
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
    trace1s=stream1s[t]
    trace2s=stream2s[t]
    trace3s=stream3s[t]
    trace4s=stream4s[t]
    trace5s=stream5s[t]
    trace6s=stream6s[t]
 
    st1=trace1.data
    st2=trace2.data
    st3=trace3.data
    st4=trace4.data
    st5=trace5.data
    st6=trace6.data
    st1s=trace1s.data
    st2s=trace2s.data
    st3s=trace3s.data
    st4s=trace4s.data
    st5s=trace5s.data
    st6s=trace6s.data
    
    cft1=recursive_sta_lta(st1, nsta, nlta)
    cft2=recursive_sta_lta(st2, nsta, nlta)
    cft3=recursive_sta_lta(st3, nsta, nlta)
    cft4=recursive_sta_lta(st4, nsta, nlta)
    cft5=recursive_sta_lta(st5, nsta, nlta)
    cft6=recursive_sta_lta(st6, nsta, nlta)
    cft1s=recursive_sta_lta(st1s, nsta, nlta)
    cft2s=recursive_sta_lta(st2s, nsta, nlta)
    cft3s=recursive_sta_lta(st3s, nsta, nlta)
    cft4s=recursive_sta_lta(st4s, nsta, nlta)
    cft5s=recursive_sta_lta(st5s, nsta, nlta)
    cft6s=recursive_sta_lta(st6s, nsta, nlta)
                                   
#    plot_trigger(trs2, cft2, trig_on, trig_off) 
    on_off_all=np.zeros(shape=(1,14))
    
    on_off1 = trigger_onset(cft1,trig_on,trig_off)
    on_off2 = trigger_onset(cft2,trig_on,trig_off)
    on_off3 = trigger_onset(cft3,trig_on,trig_off)
    on_off4 = trigger_onset(cft4,trig_on,trig_off)
    on_off5 = trigger_onset(cft5,trig_on,trig_off)
    on_off6 = trigger_onset(cft6,trig_on,trig_off)
    on_off1s = trigger_onset(cft1s,trig_on,trig_off)
    on_off2s = trigger_onset(cft2s,trig_on,trig_off)
    on_off3s = trigger_onset(cft3s,trig_on,trig_off)
    on_off4s = trigger_onset(cft4s,trig_on,trig_off)
    on_off5s = trigger_onset(cft5s,trig_on,trig_off)
    on_off6s = trigger_onset(cft6s,trig_on,trig_off)
    
#%% All LB01 triggered events
    event_day=[]

    for x in range(0,len(on_off1)):
        start=trace1.stats.starttime
        trace1.detrend(type='linear')
        trace1.detrend(type='demean')
        trace1.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset1=start+(on_off1[x,0]/sr) 
        offset1= start+(on_off1[x,1]/sr) 
        event_len = offset1-onset1
        tr = trace1.slice(starttime=onset1 - 15 , endtime=onset1 + 55)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        
        rt=onset1.timestamp
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp-rt) > 40:
        
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
             
            if (4 < cf < 11) and (3 < peak < 12) and (2 < bwid < 12) and (4 < bwid25) and(event_len > 18) :
            
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
        trace2.detrend(type='linear')
        trace2.detrend(type='demean')
        trace2.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset2=start2+(on_off2[x,0]/sr) 
        offset2= start2+(on_off2[x,1]/sr)
        event_len = offset2-onset2
        tr = trace2.slice(starttime=onset2 -15 , endtime=onset2 + 55)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        
        rt=onset2.timestamp	
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp-rt) > 40: 
        
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
    
            if (4 < cf < 10) and (3 < peak < 13) and (1.5 < bwid < 10) and (4 < bwid25 )  and (event_len > 18) :
                
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
        trace3.detrend(type='linear')
        trace3.detrend(type='demean')
        trace3.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset3=start3+(on_off3[x,0]/sr) 
        offset3= start3+(on_off3[x,1]/sr)
        event_len = offset3-onset3
        tr = trace3.slice(starttime=onset3 -15, endtime=onset3 + 55 )
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime
        rt=onset3.timestamp	
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp-rt) > 40:
 
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
        
            if (4 < cf < 11) and (3 < peak < 12) and (1.2 < bwid < 9) and (4 < bwid25 )  and (event_len > 18) :
                
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
        trace4.detrend(type='linear')
        trace4.detrend(type='demean')
        trace4.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset4=start4+(on_off4[x,0]/sr) 
        offset4= start4+(on_off4[x,1]/sr) 
        event_len = offset4-onset4
        tr = trace4.slice(starttime=onset4 -15, endtime=onset4 + 55)        
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        
        rt=onset4.timestamp	
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp-rt) > 40: 
            
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
        
            if (4 < cf < 10) and (2.25 < peak < 10) and (2 < bwid < 10) and (4 < bwid25 )  and (event_len > 18) :
                           
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
        trace5.detrend(type='linear')
        trace5.detrend(type='demean')
        trace5.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset5=start5+(on_off5[x,0]/sr) 
        offset5= start5+(on_off5[x,1]/sr) 
        event_len = offset5-onset5
        tr = trace5.slice(starttime=onset5 -15, endtime=onset5 + 55)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime
        
        rt=onset5.timestamp	
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp-rt) > 40:
                
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
     
            if (4.5 < cf < 9) and (3 < peak < 9) and (2 < bwid < 12) and (2.5 < bwid25 )  and (event_len > 18) :
            
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
        trace6.detrend(type='linear')
        trace6.detrend(type='demean')
        trace6.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset6=start6+(on_off6[x,0]/sr) 
        offset6= start6+(on_off6[x,1]/sr)
        event_len = offset6-onset6
        tr = trace6.slice(starttime=onset6 -15, endtime=onset6 + 55)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        
        rt=onset6.timestamp	
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp-rt) > 40:
                
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
        
            if (5 < cf < 9) and (3 < peak < 9) and (2.5 < bwid < 9) and (4 < bwid25 )  and (event_len > 18):
            
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

#%% All LS01 triggered events

    for x in range(0,len(on_off1s)):
        start=trace1s.stats.starttime
        trace1s.detrend(type='linear')
        trace1s.detrend(type='demean')
        trace1s.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset1s=start+(on_off1s[x,0]/sr) 
        offset1s= start+(on_off1s[x,1]/sr) 
        event_len = offset1s-onset1s
        tr = trace1s.slice(starttime=onset1s - 15 , endtime=onset1s + 55)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        
        rt=onset1s.timestamp
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp-rt) > 40:
                
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
            
            if (3.5 < cf < 12) and (3 < peak < 9) and ( bwid < 12)  and (3.7 < bwid25) and (event_len > 18) :
            
                trs_e = obspy.signal.filter.envelope(tr.data)
                top_v_eq,top_eq,corell_eq = corel(trc2_s1,trs_e,shift) 
                 
                if abs(top_v_eq) > crite:
                
                    file_ls01.write(str(onset1s))
                    file_ls01.write("\n")                                            
                    event_day.append(onset1s)
                	
                    jd=onset1s.julday
                    yr=onset1s.year
                    mo=onset1s.month
                    da=onset1s.day
                    hr=onset1s.hour
                    mi=onset1s.minute
                    se=onset1s.second
                    ms=onset1s.microsecond/10000

                    value = rt
                    near,ind=find_nearest(on_off_all[:,0], value)
    
                    if abs(value-near) > 30: #if outside 30s of an event, save to full list  
                        on_off_all[event_count][0]=rt
                        on_off_all[event_count][1]=event_len
                        on_off_all[event_count][9]=1
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                        event_day.append(onset1s)
                    else: 
                        on_off_all[ind][9]=1
                        
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
                                    corr_event2[corr_event_count2][14]=on_off_all[ind][7]
                                    corr_event2[corr_event_count2][15]=on_off_all[ind][8]
                                    corr_event_count2 += 1
                                    corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    corr_day2.append(onset1s) 
                            else:
                                corr_event2[indx][15]=1                        
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
                            corr_event2[corr_event_count2][15]=1   
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset1s)
#                            
    #%% LS02 only events

    for x in range(0,len(on_off2s)):
        start2s=trace2s.stats.starttime
        trace2s.detrend(type='linear')
        trace2s.detrend(type='demean')
        trace2s.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset2s=start2s+(on_off2s[x,0]/sr) 
        offset2s= start2s+(on_off2s[x,1]/sr)
        event_len = offset2s-onset2s
        tr = trace2s.slice(starttime=onset2s -15 , endtime=onset2s + 55)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        
        rt=onset2s.timestamp	
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp-rt) > 40: 
                
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
    
            if (6 < cf < 13) and (3 < peak < 15) and ( bwid < 13) and (6 < bwid25 )  and (event_len > 18) :
                   
                trs_e = obspy.signal.filter.envelope(tr.data)
                top_v_eq,top_eq,corell_eq = corel(trc2_s2,trs_e,shift)
    
                if abs(top_v_eq) > crite:
    
                    file_ls02.write(str(onset2s))
                    file_ls02.write("\n")
                    
                    jd=onset2s.julday
                    yr=onset2s.year
                    mo=onset2s.month
                    da=onset2s.day
                    hr=onset2s.hour
                    mi=onset2s.minute
                    se=onset2s.second
                    ms=onset2s.microsecond/10000                                            
    
                    value = rt
                    near,ind=find_nearest(on_off_all[:,0], value)
    
                    if abs(value-near) > 30: #if outside 30s of an event, save to full list  
                        on_off_all[event_count][0]=rt
                        on_off_all[event_count][1]=event_len
                        on_off_all[event_count][9]=1
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                        event_day.append(onset2s)
                    else: 
                        on_off_all[ind][9]=1
                        
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
                                    corr_event2[corr_event_count2][14]=on_off_all[ind][7]
                                    corr_event2[corr_event_count2][15]=on_off_all[ind][8]
                                    corr_event2[corr_event_count2][16]=on_off_all[ind][9]
                                    corr_event_count2 += 1
                                    corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    corr_day2.append(onset2s) 
                            else:
                                corr_event2[indx][16]=1                        
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
                            corr_event2[corr_event_count2][16]=1   
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset2s)

#%% LS03 only events

    for x in range(0,len(on_off3s)):
        start3s=trace3s.stats.starttime
        trace3s.detrend(type='linear')
        trace3s.detrend(type='demean')
        trace3s.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset3s=start3s+(on_off3s[x,0]/sr) 
        offset3s= start3s+(on_off3s[x,1]/sr)
        event_len = offset3s-onset3s
        tr = trace3s.slice(starttime=onset3s -15, endtime=onset3s + 55 )
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime
        
        rt=onset3s.timestamp	
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp-rt) > 40:
                
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
        
            if (4.5 < cf < 10.5) and (3 < peak < 10) and (2 < bwid < 15) and (3 < bwid25 )  and (event_len > 18)  :
            
                trs_e = obspy.signal.filter.envelope(tr.data)          
                top_v_eq,top_eq,corell_eq = corel(trc2_s3,trs_e,shift)            
           
                if abs(top_v_eq) > crite:
                    
                    file_ls03.write(str(onset3s))
                    file_ls03.write("\n") 
                    
                    jd=onset3s.julday
                    yr=onset3s.year
                    mo=onset3s.month
                    da=onset3s.day
                    hr=onset3s.hour
                    mi=onset3s.minute
                    se=onset3s.second
                    ms=onset3s.microsecond/10000
                    
                    value = rt
                    near,ind=find_nearest(on_off_all[:,0], value)
    
                    if abs(value-near) > 30: #if outside 30s of an event, save to full list  
                        on_off_all[event_count][0]=rt
                        on_off_all[event_count][1]=event_len
                        on_off_all[event_count][10]=1
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                        event_day.append(onset3s)
                    else: 
                        on_off_all[ind][10]=1
                        
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
                                    corr_event2[corr_event_count2][14]=on_off_all[ind][7]
                                    corr_event2[corr_event_count2][15]=on_off_all[ind][8]
                                    corr_event2[corr_event_count2][16]=on_off_all[ind][9]
                                    corr_event2[corr_event_count2][17]=on_off_all[ind][10]
                                    corr_event_count2 += 1
                                    corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    corr_day2.append(onset3s) 
                            else:
                                corr_event2[indx][17]=1                        
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
                            corr_event2[corr_event_count2][17]=1   
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset3s) 
                            
    #%% LS04 only events

    for x in range(0,len(on_off4s)):
        start4s=trace4s.stats.starttime
        trace4s.detrend(type='linear')
        trace4s.detrend(type='demean')
        trace4s.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset4s=start4s+(on_off4s[x,0]/sr) 
        offset4s= start4s+(on_off4s[x,1]/sr) 
        event_len = offset4s-onset4s
        tr = trace4s.slice(starttime=onset4s -15, endtime=onset4s + 55)        
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime
        
        rt=onset4s.timestamp	
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp-rt) > 40:   
                
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
        
            if (4 < cf < 10) and (3 < peak < 10) and ( 2 < bwid < 12) and (4 < bwid25 )  and (event_len > 18) :
                         
                trs_e = obspy.signal.filter.envelope(tr.data)
                top_v_eq,top_eq,corell_eq = corel(trc2_s4,trs_e,shift)
                               
                if abs(top_v_eq) > crite:
                            
                    file_ls04.write(str(onset4s))
                    file_ls04.write("\n") 
                    	
                    jd=onset4s.julday
                    yr=onset4s.year
                    mo=onset4s.month
                    da=onset4s.day
                    hr=onset4s.hour
                    mi=onset4s.minute
                    se=onset4s.second
                    ms=onset4s.microsecond/10000
    
                    value = rt
                    near,ind=find_nearest(on_off_all[:,0], value)
    #                    corr_near,indx=find_nearest(corr_event2[:,0], value)
    
                    if abs(value-near) > 30: #if outside 30s of an event, save to full list  
                        on_off_all[event_count][0]=rt
                        on_off_all[event_count][1]=event_len
                        on_off_all[event_count][11]=1
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                        event_day.append(onset4s)
                    else: 
                        on_off_all[ind][11]=1
                        
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
                                    corr_event2[corr_event_count2][14]=on_off_all[ind][7]
                                    corr_event2[corr_event_count2][15]=on_off_all[ind][8]
                                    corr_event2[corr_event_count2][16]=on_off_all[ind][9]
                                    corr_event2[corr_event_count2][17]=on_off_all[ind][10]
                                    corr_event2[corr_event_count2][18]=on_off_all[ind][11]
                                    corr_event_count2 += 1
                                    corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    corr_day2.append(onset4s) 
                            else:
                                corr_event2[indx][18]=1                        
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
                            corr_event2[corr_event_count2][18]=1   
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset4s)

    #%% LS05 only events

    for x in range(0,len(on_off5s)):
        start5s=trace5s.stats.starttime
        trace5s.detrend(type='linear')
        trace5s.detrend(type='demean')
        trace5s.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset5s=start5s+(on_off5s[x,0]/sr) 
        offset5s= start5s+(on_off5s[x,1]/sr) 
        event_len = offset5s-onset5s
        tr = trace5s.slice(starttime=onset5s -15, endtime=onset5s + 55)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        
        rt=onset5s.timestamp	
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp-rt) > 40:
                
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
     
            if (5 < cf < 11) and (5 < peak < 12) and (3 < bwid < 10) and (5 < bwid25 )  and (event_len > 18) :# 
            
                trs_e = obspy.signal.filter.envelope(tr.data)
                top_v_eq,top_eq,corell_eq = corel(trc2_s5,trs_e,shift)
           
                if abs(top_v_eq) > crite:
                    
                    file_ls05.write(str(onset5s))
                    file_ls05.write("\n") 
                    
                    jd=onset5s.julday
                    yr=onset5s.year
                    mo=onset5s.month
                    da=onset5s.day
                    hr=onset5s.hour
                    mi=onset5s.minute
                    se=onset5s.second
                    ms=onset5s.microsecond/10000
                    
                    value = rt
                    near,ind=find_nearest(on_off_all[:,0], value)
    #                    corr_near,indx=find_nearest(corr_event2[:,0], value)
    
                    if abs(value-near) > 30: #if outside 30s of an event, save to full list  
                        on_off_all[event_count][0]=rt
                        on_off_all[event_count][1]=event_len
                        on_off_all[event_count][12]=1
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                        event_day.append(onset5s)
                    else: 
                        on_off_all[ind][12]=1
                        
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
                                    corr_event2[corr_event_count2][14]=on_off_all[ind][7]
                                    corr_event2[corr_event_count2][15]=on_off_all[ind][8]
                                    corr_event2[corr_event_count2][16]=on_off_all[ind][9]
                                    corr_event2[corr_event_count2][17]=on_off_all[ind][10]
                                    corr_event2[corr_event_count2][18]=on_off_all[ind][11]
                                    corr_event2[corr_event_count2][19]=on_off_all[ind][12]
                                    corr_event_count2 += 1
                                    corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    corr_day2.append(onset5s)
                            else:
                                corr_event2[indx][19]=1
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
                            corr_event2[corr_event_count2][19]=1   
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset5s)
          
    #%% LS06 only events
                        
    for x in range(0,len(on_off6s)):
        start6s=trace6s.stats.starttime
        trace6s.detrend(type='linear')
        trace6s.detrend(type='demean')
        trace6s.filter("bandpass", freqmin=fmin,freqmax=fmax)
        onset6s=start6s+(on_off6s[x,0]/sr) 
        offset6s= start6s+(on_off6s[x,1]/sr)
        event_len = offset6s-onset6s
        tr = trace6s.slice(starttime=onset6s -15, endtime=onset6s + 55)
        tr.detrend(type='demean')
        tr.detrend(type='linear')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        
        rt=onset6s.timestamp	
        near_exp,ix=find_nearest(ts[:,0], rt)
        if abs(near_exp - rt) > 40:
                
            peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
        
            if (5 < cf < 11) and (3.5 < peak < 10) and (1 < bwid < 10) and (5 < bwid25 )  and (event_len > 18) :
            
                trs_e = obspy.signal.filter.envelope(tr.data)
                top_v_eq,top_eq,corell_eq = corel(trc2_s6,trs_e,shift)
                            
                if abs(top_v_eq) > crite:
                 
                    file_ls06.write(str(onset6s))
                    file_ls06.write("\n") 
                    
                    jd=onset6s.julday
                    yr=onset6s.year
                    mo=onset6s.month
                    da=onset6s.day
                    hr=onset6s.hour
                    mi=onset6s.minute
                    se=onset6s.second
                    ms=onset6s.microsecond/10000
                    
                    value = rt
                    near,ind=find_nearest(on_off_all[:,0], value)
    #                    corr_near,indx=find_nearest(corr_event2[:,0], value)
    
                    if abs(value-near) > 30: #if outside 30s of an event, save to full list  
                        on_off_all[event_count][0]=rt
                        on_off_all[event_count][1]=event_len
                        on_off_all[event_count][13]=1
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                        event_day.append(onset6s)
                    else: 
                        on_off_all[ind][13]=1
                        
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
                                    corr_event2[corr_event_count2][15]=on_off_all[ind][8]
                                    corr_event2[corr_event_count2][16]=on_off_all[ind][9]
                                    corr_event2[corr_event_count2][17]=on_off_all[ind][10]
                                    corr_event2[corr_event_count2][18]=on_off_all[ind][11]
                                    corr_event2[corr_event_count2][19]=on_off_all[ind][12]
                                    corr_event2[corr_event_count2][20]=on_off_all[ind][13]
                                    corr_event_count2 += 1
                                    corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                    corr_day2.append(onset6s)
                            else:
                                corr_event2[indx][20]=1
                    
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
                            corr_event2[corr_event_count2][20]=1   
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset6s)                             
                                
#%%            
    event_day.sort()
    corr_day2.sort()  

    for p in range(0,len(corr_day2)):
        file_cor2.write(str(corr_day2[p]))
        file_cor2.write("\n")  
            

#%%   
print('Saving Data')          
col=0       
corr_event2=corr_event2[np.argsort(corr_event2[:,col])] 
if corr_event2[0][0]== 0.0:
    corr_event2 = np.delete(corr_event2, (0), axis=0)       
np.savetxt("/Users/william/Documents/scanner/all_stations/EQ_all_coincidence_2_month_30.csv", corr_event2,delimiter=",",header="Time_stamp,Event_length,Day_number,Year,Day_of_year,Hour,Min,Sec,Milisec,LBO1,LB02,LB03,LB04,LB05,LB06,LSO1,LS02,LS03,LS04,LS05,LS06")

data_stream = Stream()
for r in range(0,len(corr_event2)):
    if corr_event2[r][9]==1:
        tr1a=stream1[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr1_save = tr1a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr1_save.detrend(type='demean')
        tr1_save.detrend(type='linear')
        if type(tr1_save[0]) == np.float64:
            data_stream.append(tr1_save)
    if corr_event2[r][10]==1:
        tr2a=stream2[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr2_save = tr2a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr2_save.detrend(type='demean')
        tr2_save.detrend(type='linear')
        if type(tr2_save[0]) == np.float64:
            data_stream.append(tr2_save)    
    if corr_event2[r][11]==1:
        tr3a=stream3[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr3_save = tr3a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr3_save.detrend(type='demean')
        tr3_save.detrend(type='linear')
        if type(tr3_save[0]) == np.float64:
            data_stream.append(tr3_save)         
    if corr_event2[r][12]==1:
        tr4a=stream4[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr4_save = tr4a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr4_save.detrend(type='demean')
        tr4_save.detrend(type='linear')
        if type(tr4_save[0]) == np.float64:
            data_stream.append(tr4_save)    
    if corr_event2[r][13]==1:
        tr5a=stream5[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr5_save = tr5a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr5_save.detrend(type='demean')
        tr5_save.detrend(type='linear')
        if type(tr5_save[0]) == np.float64:
            data_stream.append(tr5_save)
    if corr_event2[r][14]==1:
        tr6a=stream6[int(corr_event2[r][2])] 
        dstart=UTCDateTime(corr_event2[r][0])
        tr6_save = tr6a.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr6_save.detrend(type='demean')
        tr6_save.detrend(type='linear')
        if type(tr6_save[0]) == np.float64:
            data_stream.append(tr6_save)
    if corr_event2[r][15]==1:
        tr1sa=stream1s[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr1s_save = tr1sa.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr1s_save.detrend(type='demean')
        tr1s_save.detrend(type='linear')
        if type(tr1s_save[0]) == np.float64:
            data_stream.append(tr1s_save)
    if corr_event2[r][16]==1:
        tr2sa=stream2s[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr2s_save = tr2sa.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr2s_save.detrend(type='demean')
        tr2s_save.detrend(type='linear')
        if tr2s_save.stats.npts > 1:
            if type(tr2s_save[0]) == np.float64 :
                data_stream.append(tr2s_save)    
    if corr_event2[r][17]==1:
        tr3sa=stream3s[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr3s_save = tr3sa.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr3s_save.detrend(type='demean')
        tr3s_save.detrend(type='linear')
        if type(tr3s_save[0]) == np.float64:
            data_stream.append(tr3s_save)         
    if corr_event2[r][18]==1:
        tr4sa=stream4s[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr4s_save = tr4sa.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr4s_save.detrend(type='demean')
        tr4s_save.detrend(type='linear')
        if type(tr4s_save[0]) == np.float64:
            data_stream.append(tr4s_save)    
    if corr_event2[r][19]==1:
        tr5sa=stream5s[int(corr_event2[r][2])]
        dstart=UTCDateTime(corr_event2[r][0])
        tr5s_save = tr5sa.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr5s_save.detrend(type='demean')
        tr5s_save.detrend(type='linear')
        if type(tr5s_save[0]) == np.float64:
            data_stream.append(tr5s_save)
    if corr_event2[r][20]==1:
        tr6sa=stream6s[int(corr_event2[r][2])] 
        dstart=UTCDateTime(corr_event2[r][0])
        tr6s_save = tr6sa.slice(starttime=dstart - 10, endtime=dstart+corr_event2[r][1]+20)
        tr6s_save.detrend(type='demean')
        tr6s_save.detrend(type='linear')
        if type(tr6s_save[0]) == np.float64:
            data_stream.append(tr6s_save)

if len(data_stream) > 0:   
    data_stream.write('/Users/william/Documents/scanner/output_data/EQ_all_data_stream_2_month_30.mseed', format='MSEED')         

        
#file_sta.close()
file_cor2 .close()
file_lb01.close()
file_lb02.close()
file_lb03.close()
file_lb04.close()
file_lb05.close()
file_lb06.close()
file_ls01.close()
file_ls02.close()
file_ls03.close()
file_ls04.close()
file_ls05.close()
file_ls06.close()

print('end of scan')
T2=time.clock()
elapsed= T2-T1
print('Time taken:', elapsed)       


