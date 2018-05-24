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
   #%% constants  
shift=500               # was 200 but increased to 500 to get higher Xcorr - will be slower
crite=0.6               # V8 =0.55 V8b = 0.6   
window=30
fmin=3
fmax=12
#%%
file_sta = open('/Users/william/Documents/scanner/all_stations/full_eq_scan_v8b.txt','w')
file_cor2 = open('/Users/william/Documents/scanner/all_stations/full_correlated2_eq_scan_v8b.txt','w')
file_cor3 = open('/Users/william/Documents/scanner/all_stations/full_correlated3_eq_scan_v8b.txt','w')
file_lb01= open('/Users/william/Documents/scanner/LB01/LB01_full_eq_scan_v8b.txt','w')
file_lb02= open('/Users/william/Documents/scanner/LB02/LB02_full_eq_scan_v8b.txt','w')
file_lb03= open('/Users/william/Documents/scanner/LB03/LB03_full_eq_scan_v8b.txt','w')
file_lb04= open('/Users/william/Documents/scanner/LB04/LB04_full_eq_scan_v8b.txt','w')
file_lb05= open('/Users/william/Documents/scanner/LB05/LB05_full_eq_scan_v8b.txt','w')
file_lb06= open('/Users/william/Documents/scanner/LB06/LB06_full_eq_scan_v8b.txt','w')
file_lb07= open('/Users/william/Documents/scanner/LB07/LB07_full_eq_scan_v8b.txt','w')

#%% Read in waveforms 
  
trc2_e1,trc2_e2,trc2_e3,trc2_e4,trc2_e5,trc2_e6,trc2_e7=stacked_LBz_eq() 

#%% open empty streams for each station

stream1 = Stream()
stream2 = Stream()
stream3 = Stream()
stream4 = Stream()
stream5 = Stream()
stream6 = Stream()
stream7 = Stream()

#%% Read in Each station

num_active=[]

for x in range(0,10):
    
    seis1,seis2,seis3,seis4,seis5,seis6,seis7,num=get_LBz(x)
    
    num_active.append(num)

    tr1=seis1[0]
    tr2=seis2[0]
    tr3=seis3[0]
    tr4=seis4[0]
    tr5=seis5[0]
    tr6=seis6[0]
    tr7=seis7[0]

    stream1.append(tr1)
    stream2.append(tr2)
    stream3.append(tr3)
    stream4.append(tr4)
    stream5.append(tr5)
    stream6.append(tr6)
    stream7.append(tr7)
   
print('files read in') 

#%%
sr = seis1[0].stats.sampling_rate
nsta=int(1*sr)                                      #2
nlta=int(20*sr)                                     #20
trig_on=7                                            #8
trig_off=0.2                                         #0.2

count=0

for t in range(0,len(stream1),1): 
    event_count=0
    corr_event_count2=0
    corr_event_count3=0
    print('Day',t+1,'of',len(stream1))
    file_sta.write("\n")
    file_cor2.write("\n")
    file_cor3.write("\n")
    corr_event2=np.zeros(shape=(1,2))
    corr_day2=[]
    corr_event3=np.zeros(shape=(1,2))
    corr_day3=[]
    
    trace1=stream1[t]
    trace2=stream2[t]
    trace3=stream3[t]
    trace4=stream4[t]
    trace5=stream5[t]
    trace6=stream6[t]
    trace7=stream7[t]
 
    st1=trace1.data
    st2=trace2.data
    st3=trace3.data
    st4=trace4.data
    st5=trace5.data
    st6=trace6.data
    st7=trace7.data
    
    cft1=recursive_sta_lta(st1, nsta, nlta)
    cft2=recursive_sta_lta(st2, nsta, nlta)
    cft3=recursive_sta_lta(st3, nsta, nlta)
    cft4=recursive_sta_lta(st4, nsta, nlta)
    cft5=recursive_sta_lta(st5, nsta, nlta)
    cft6=recursive_sta_lta(st6, nsta, nlta)
    cft7=recursive_sta_lta(st7, nsta, nlta)
                                   
#    plot_trigger(trs2, cft2, trig_on, trig_off) 
    on_off_all=np.zeros(shape=(1,2))
    
    on_off1 = trigger_onset(cft1,trig_on,trig_off)
    on_off2 = trigger_onset(cft2,trig_on,trig_off)
    on_off3 = trigger_onset(cft3,trig_on,trig_off)
    on_off4 = trigger_onset(cft4,trig_on,trig_off)
    on_off5 = trigger_onset(cft5,trig_on,trig_off)
    on_off6 = trigger_onset(cft6,trig_on,trig_off)
    on_off7 = trigger_onset(cft7,trig_on,trig_off)

#    t0 = UTCDateTime(2014, 11, 24, 0, 0, 0) # start of first day read in
#    start = t0 + t*24*60*60

#%% All LB01 triggered events
    event_day=[]

    for x in range(0,len(on_off1)):
        start=trace1.stats.starttime
        tr = trace1.slice(starttime=start+(on_off1[x,0]/sr)-15 , endtime=start+(on_off1[x,1]/sr)+25)
        tr_len=tr.stats.endtime - tr.stats.starttime 
        onset1=start+(on_off1[x,0]/sr)
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
        
        tr.detrend(type='demean')
        amp=abs(tr.max())
        whole = sum(abs(tr.data))
        av_amp = whole/(100*tr_len)
        amp_rat = amp/av_amp
        # only look for things above noise level 
        
#        if amp > 300:
#            if amp_rat < 15:
        if 5 < cf < 11 :  
            if 3 < peak < 13:
                if 2 < bwid < 12:
                    
                    trs_e = obspy.signal.filter.envelope(tr.data)
    
                    top_v_eq,top_eq,corell_eq = corel(trc2_e1,trs_e,shift)                    
                     
                    if abs(top_v_eq) > crite:

                        file_lb01.write(str(onset1))
                        file_lb01.write("\n")                                            
                        event_day.append(onset1)
                        
                        on_off_all[event_count][0]=on_off1[x,0]
                        on_off_all[event_count][1]=on_off1[x,1]
                        event_count += 1
                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))     

    #%% LB02 only events

    for x in range(0,len(on_off2)):
        start=trace2.stats.starttime
        tr = trace2.slice(starttime=start+(on_off2[x,0]/sr)-15 , endtime=start+(on_off2[x,1]/sr)+25)
        tr.detrend(type='demean')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        onset2=start+(on_off2[x,0]/sr)
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      

    
        amp=abs(tr.max())
        whole = sum(abs(tr.data))
        av_amp = whole/(100*tr_len)
        amp_rat = amp/av_amp
        # only look for things above noise level 
        
#        if amp > 3000:
#            if amp_rat < 15:
        if 5 < cf < 11 :  
            if 4 < peak < 13:
                if 0.2 < bwid < 9:
#                            
                    trs_e = obspy.signal.filter.envelope(tr.data)
                  
                    top_v_eq,top_eq,corell_eq = corel(trc2_e2,trs_e,shift)

                    if abs(top_v_eq) > crite:

                        file_lb02.write(str(onset2))
                        file_lb02.write("\n")                                            

                        value = on_off2[x,0] 
                        near,ind=find_nearest(on_off_all[:,0], value)
                
                        if abs(value-near) > 100*30: #if outside 30s of an event, save to full list  
                            on_off_all[event_count][0]=on_off2[x][0]
                            on_off_all[event_count][1]=on_off2[x][1]
                            event_count += 1
                            on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                            event_day.append(onset2)
                        else: # if already found, save as a 2x coincidence event
                            corr_event2[corr_event_count2][0]=on_off2[x][0]
                            corr_event2[corr_event_count2][1]=on_off2[x][1]
                            corr_event_count2 += 1
                            corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day2.append(onset2)               

#%% LB03 only events

    for x in range(0,len(on_off3)):
        start=trace3.stats.starttime
        tr = trace3.slice(starttime=start+(on_off3[x,0]/sr)-15 , endtime=start+(on_off3[x,1]/sr)+25)
        tr.detrend(type='demean')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        onset3=start+(on_off3[x,0]/sr)
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
    
        amp=abs(tr.max())
        whole = sum(abs(tr.data))
        av_amp = whole/(100*tr_len)
        amp_rat = amp/av_amp
        # only look for things above noise level 
        
#        if amp > 40000:
#            if amp_rat < 15:
        if 5 < cf < 12.5 :  
            if 4 < peak < 13.5:
                if 1 < bwid < 10.5:
  
                    trs_e = obspy.signal.filter.envelope(tr.data)
                  
                    top_v_eq,top_eq,corell_eq = corel(trc2_e3,trs_e,shift)
                    
               
                    if abs(top_v_eq) > crite:
                        
                        file_lb03.write(str(onset3))
                        file_lb03.write("\n") 
                        
                        value = on_off3[x,0] 
                        near,ind=find_nearest(on_off_all[:,0], value)
                        corr_near,indx=find_nearest(corr_event2[:,0], value)
                        
                        if abs(value-corr_near) < 100*30: # if within 30s of value in 2x coincidence - now in 3x
                            corr_event3[corr_event_count3][0]=on_off3[x][0]
                            corr_event3[corr_event_count3][1]=on_off3[x][1]
                            corr_event_count3 += 1
                            corr_event3 = np.lib.pad(corr_event3, ((0,1),(0,0)), 'constant', constant_values=(0))
                            corr_day3.append(onset3)
                        else: #  not within 30s of 2x
                            if abs(value-near) > 100*30: # not in list of all events - add to it 
                                on_off_all[event_count][0]=on_off3[x][0]
                                on_off_all[event_count][1]=on_off3[x][1]
                                event_count += 1
                                on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                                event_day.append(onset3)
                            else: # add to 2x coincidence if it is in events list only
                                corr_event2[corr_event_count2][0]=on_off3[x][0]
                                corr_event2[corr_event_count2][1]=on_off3[x][1]
                                corr_event_count2 += 1
                                corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                corr_day2.append(onset3)
                                        
    #%% LB04 only events

    for x in range(0,len(on_off4)):
        start=trace4.stats.starttime
        tr = trace4.slice(starttime=start+(on_off4[x,0]/sr)-15 , endtime=start+(on_off4[x,1]/sr)+25)        
        tr.detrend(type='demean')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        onset4=start+(on_off4[x,0]/sr)
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
    
        amp=abs(tr.max())
        whole = sum(abs(tr.data))
        av_amp = whole/(100*tr_len)
        amp_rat = amp/av_amp
        # only look for things above noise level 
        
#        if amp > 20000:
#            if amp_rat < 15:
        if 3.5 < cf < 11 :  
            if 2 < peak < 12:
                if 0.7 < bwid < 11:
#                           
                    trs_e = obspy.signal.filter.envelope(tr.data)
                  
                    top_v_eq,top_eq,corell_eq = corel(trc2_e4,trs_e,shift)
                                   
                    if abs(top_v_eq) > crite:
                                
                        file_lb04.write(str(onset4))
                        file_lb04.write("\n") 
                        
                        value = on_off4[x,0] 
                        near,ind=find_nearest(on_off_all[:,0], value)
                        corr_near,indx=find_nearest(corr_event2[:,0], value)
                        corr_near3,indx3=find_nearest(corr_event3[:,0], value)
                        
                        if abs(value-corr_near) < 100*30:
                            if abs(value-corr_near3) > 30*100: # add to 3x if in 2x not already in 3x
                                corr_event3[corr_event_count3][0]=on_off4[x][0]
                                corr_event3[corr_event_count3][1]=on_off4[x][1]
                                corr_event_count3 += 1
                                corr_event3 = np.lib.pad(corr_event3, ((0,1),(0,0)), 'constant', constant_values=(0))
                                corr_day3.append(onset4)
                        else: 
                            if abs(value-near) > 100*30: #15s  
                                on_off_all[event_count][0]=on_off4[x][0]
                                on_off_all[event_count][1]=on_off4[x][1]
                                event_count += 1
                                on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                                event_day.append(onset4)
                            else:
                                corr_event2[corr_event_count2][0]=on_off4[x][0]
                                corr_event2[corr_event_count2][1]=on_off4[x][1]
                                corr_event_count2 += 1
                                corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                corr_day2.append(onset4)

    #%% LB05 only events

    for x in range(0,len(on_off5)):
        start=trace5.stats.starttime
        tr = trace5.slice(starttime=start+(on_off5[x,0]/sr)-15 , endtime=start+(on_off5[x,1]/sr)+25)
        tr.detrend(type='demean')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        onset5=start+(on_off5[x,0]/sr)
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
    
        amp=abs(tr.max())
        whole = sum(abs(tr.data))
        av_amp = whole/(100*tr_len)
        amp_rat = amp/av_amp
        # only look for things above noise level 
        
#        if amp > 20000:
#            if amp_rat < 15:
        if 4 < cf < 10 :  
            if 2.5 < peak < 11:
                if 1.5 < bwid < 10.5:
                 
                    trs_e = obspy.signal.filter.envelope(tr.data)
                  
                    top_v_eq,top_eq,corell_eq = corel(trc2_e5,trs_e,shift)
               
                    if abs(top_v_eq) > crite:
                        
                        file_lb05.write(str(onset5))
                        file_lb05.write("\n") 
                        
                        value = on_off5[x,0] 
                        near,ind=find_nearest(on_off_all[:,0], value)
                        corr_near,indx=find_nearest(corr_event2[:,0], value)
                        corr_near3,indx3=find_nearest(corr_event3[:,0], value)
                        
                        if abs(value-corr_near) < 100*30:
                            if abs(value-corr_near3) > 30*100: # add to 3x if in 2x not already in 3x                        
                                corr_event3[corr_event_count3][0]=on_off5[x][0]
                                corr_event3[corr_event_count3][1]=on_off5[x][1]
                                corr_event_count3 += 1
                                corr_event3 = np.lib.pad(corr_event3, ((0,1),(0,0)), 'constant', constant_values=(0))
                                corr_day3.append(onset5)
                        else:
                            if abs(value-near) > 100*30: #15s  
                                on_off_all[event_count][0]=on_off5[x][0]
                                on_off_all[event_count][1]=on_off5[x][1]
                                event_count += 1
                                on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                                event_day.append(onset5)
                            else:
                                corr_event2[corr_event_count2][0]=on_off5[x][0]
                                corr_event2[corr_event_count2][1]=on_off5[x][1]
                                corr_event_count2 += 1
                                corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                corr_day2.append(onset5)
                                
  
    #%% LB06 only events
                        
    for x in range(0,len(on_off6)):
        start=trace6.stats.starttime
        tr = trace6.slice(starttime=start+(on_off6[x,0]/sr)-15 , endtime=start+(on_off6[x,1]/sr)+25)
        tr.detrend(type='demean')
        tr_len=tr.stats.endtime - tr.stats.starttime 
        onset6=start+(on_off6[x,0]/sr)
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
    
        amp=abs(tr.max())
        whole = sum(abs(tr.data))
        av_amp = whole/(100*tr_len)
        amp_rat = amp/av_amp
        # only look for things above noise level 
        
#        if amp > 2000:
#            if amp_rat < 15:
        if 5 < cf < 10 :  
            if 3.5 < peak < 11:
                if 2.5 < bwid < 9.5:
                   
                    trs_e = obspy.signal.filter.envelope(tr.data)
                  
                    top_v_eq,top_eq,corell_eq = corel(trc2_e6,trs_e,shift)
                                
                    if abs(top_v_eq) > crite:
                                
                                file_lb06.write(str(onset6))
                                file_lb06.write("\n") 
                                
                                value = on_off6[x,0] 
                                near,ind=find_nearest(on_off_all[:,0], value)
                                corr_near,indx=find_nearest(corr_event2[:,0], value)
                                corr_near3,indx3=find_nearest(corr_event3[:,0], value)
                        
                                if abs(value-corr_near) < 100*30:
                                    if abs(value-corr_near3) > 30*100: # add to 3x if in 2x not already in 3x                                
                                        corr_event3[corr_event_count3][0]=on_off6[x][0]
                                        corr_event3[corr_event_count3][1]=on_off6[x][1]
                                        corr_event_count3 += 1
                                        corr_event3 = np.lib.pad(corr_event3, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        corr_day3.append(onset6)
                                else:
                                    if abs(value-near) > 100*30: #15s  
                                        on_off_all[event_count][0]=on_off6[x][0]
                                        on_off_all[event_count][1]=on_off6[x][1]
                                        event_count += 1
                                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        event_day.append(onset6)
                                    else:
                                        corr_event2[corr_event_count2][0]=on_off6[x][0]
                                        corr_event2[corr_event_count2][1]=on_off6[x][1]
                                        corr_event_count2 += 1
                                        corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        corr_day2.append(onset6)
                                        
        
    
    #%% LB07 only events                        

    for x in range(0,len(on_off7)):
        start=trace7.stats.starttime
        tr = trace7.slice(starttime=start+(on_off7[x,0]/sr)-15 , endtime=start+(on_off7[x,1]/sr)+25)
        tr.detrend(type='demean')
        tr_len=tr.stats.endtime - tr.stats.starttime
        onset7=start+(on_off7[x,0]/sr)
        
        peak,cf,bwid,bwid25 = freq_info(tr.data,tr.stats.starttime,tr.stats.endtime)      
    
        amp=abs(tr.max())
        whole = sum(abs(tr.data))
        av_amp = whole/(100*tr_len)
        amp_rat = amp/av_amp
        # only look for things above noise level 
        
#        if amp > 1000:
#            if amp_rat < 15:
        if 5.5 < cf < 10 :  
            if 4 < peak < 11:
                if 2 < bwid < 10:
                   
                    trs_e = obspy.signal.filter.envelope(tr.data)
                  
                    top_v_eq,top_eq,corell_eq = corel(trc2_e7,trs_e,shift)
                                   
                    if abs(top_v_eq) > crite:

                        file_lb07.write(str(onset7))
                        file_lb07.write("\n") 
                        
                        value = on_off7[x,0] 
                        near,ind=find_nearest(on_off_all[:,0], value)
                        corr_near,indx=find_nearest(corr_event2[:,0], value)
                        corr_near3,indx3=find_nearest(corr_event3[:,0], value)
                        
                        if abs(value-corr_near) < 100*30:
                            if abs(value-corr_near3) > 30*100: # add to 3x if in 2x not already in 3x
                                corr_event3[corr_event_count3][0]=on_off7[x][0]
                                corr_event3[corr_event_count3][1]=on_off7[x][1]
                                corr_event_count3 += 1
                                corr_event3 = np.lib.pad(corr_event3, ((0,1),(0,0)), 'constant', constant_values=(0))
                                corr_day3.append(onset7)
                        else:
                            if abs(value-near) > 100*30: #15s  
                                on_off_all[event_count][0]=on_off7[x][0]
                                on_off_all[event_count][1]=on_off7[x][1]
                                event_count += 1
                                on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))
                                event_day.append(onset7)
                            else:
                                corr_event2[corr_event_count2][0]=on_off7[x][0]
                                corr_event2[corr_event_count2][1]=on_off7[x][1]
                                corr_event_count2 += 1
                                corr_event2 = np.lib.pad(corr_event2, ((0,1),(0,0)), 'constant', constant_values=(0))
                                corr_day2.append(onset7)
    
#%%            
    event_day.sort()
    corr_day2.sort()
    corr_day3.sort()

    for p in range(0,len(event_day)):
        file_sta.write(str(event_day[p]))
        file_sta.write("\n")
    for p in range(0,len(corr_day2)):
        file_cor2.write(str(corr_day2[p]))
        file_cor2.write("\n")    
    for p in range(0,len(corr_day3)):
        file_cor3.write(str(corr_day3[p]))
        file_cor3.write("\n")   
        
file_sta.close()
file_cor2 .close()
file_cor3.close()
file_lb01.close()
file_lb02.close()
file_lb03.close()
file_lb04.close()
file_lb05.close()
file_lb06.close()
file_lb07.close()

print('end of scan')
T2=time.clock()
elapsed= T2-T1
print('Time taken:', elapsed)       

        

    
