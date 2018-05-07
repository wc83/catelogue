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
from obspy import UTCDateTime
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
file_sta = open('/Users/william/Documents/scanner/full_eq_scan_v3.txt','w')

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
trc2_e=stacked_eq(fmin,fmax)

#noise2
sample3 = read('/Users/william/Documents/lb01/14_335z.mseed')
trace3=sample3[0]
trace3.detrend(type='demean')
trace3.filter("bandpass", freqmin=fmin,freqmax=fmax)

start3= sample3[0].stats.starttime + 20*60*60 + 57*60 + 50 #time window start - several s before Pwave arrival
end3=start3 + window
trc3 = trace3.slice(starttime = start3  , endtime= end3) #cut out sample waveform with same window length as chosen event
trc3_e = obspy.signal.filter.envelope(trc3.data)
#print('reference noise(2) waveform')
#trs3.plot(type='relative',color='b', starttime=start3 , endtime=end3)

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
for x in range(0,3):
    
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

for t in range(0,len(stream1),1): 
    event_count=0
    print('Day',t+1,'of',len(stream1))
    file_sta.write("\n")
    
    trace1=stream1[t]
    trace2=stream2[t]
    trace3=stream3[t]
    trace4=stream4[t]
    trace5=stream5[t]
    trace6=stream6[t]
    trace7=stream7[t]
    trace8=stream8[t]
 
    st1=trace1.data
    st2=trace2.data
    st3=trace3.data
    st4=trace4.data
    st5=trace5.data
    st6=trace6.data
    st7=trace7.data
    st8=trace8.data
    
    cft1=recursive_sta_lta(st1, nsta, nlta)
    cft2=recursive_sta_lta(st2, nsta, nlta)
    cft3=recursive_sta_lta(st3, nsta, nlta)
    cft4=recursive_sta_lta(st4, nsta, nlta)
    cft5=recursive_sta_lta(st5, nsta, nlta)
    cft6=recursive_sta_lta(st6, nsta, nlta)
    cft7=recursive_sta_lta(st7, nsta, nlta)
    cft8=recursive_sta_lta(st8, nsta, nlta)
                                   
#    plot_trigger(trs2, cft2, trig_on, trig_off) 
    on_off_all=np.zeros(shape=(1,2))
    
    on_off1 = trigger_onset(cft1,trig_on,trig_off)
    on_off2 = trigger_onset(cft2,trig_on,trig_off)
    on_off3 = trigger_onset(cft3,trig_on,trig_off)
    on_off4 = trigger_onset(cft4,trig_on,trig_off)
    on_off5 = trigger_onset(cft5,trig_on,trig_off)
    on_off6 = trigger_onset(cft6,trig_on,trig_off)
    on_off7 = trigger_onset(cft7,trig_on,trig_off)
    on_off8 = trigger_onset(cft8,trig_on,trig_off)
    
#        #window endpoints
#    if trace1.stats.npts > 8600000:
#        start= stream1[t].stats.starttime   
#    elif trace2.stats.npts > 8600000:
#        start= stream2[t].stats.starttime  
#    if trace3.stats.npts > 8600000:
#        start= stream3[t].stats.starttime   
#    elif trace4.stats.npts > 8600000:
#        start= stream4[t].stats.starttime  
#    if trace5.stats.npts > 8600000:
#        start= stream5[t].stats.starttime   
#    elif trace6.stats.npts > 8600000:
#        start= stream6[t].stats.starttime  
#    if trace7.stats.npts > 8600000:
#        start= stream7[t].stats.starttime   
#    else:
#        break # NO station with  full day data
        
    
    t0 = UTCDateTime(2014, 11, 24, 0, 0, 0) # start of first day read in
    start = t0 + t*24*60*60

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
                          
                            top_v_eq,top_eq,corell_eq = corel(trc2_e,trs_e,shift)
                            
                            top_v_n,top_n,corell_n = corel(trc3_e,trs_e,shift)
                       
                            if abs(top_v) < crite:
                                if abs(top_v_n) > crite_n:
                                    if abs(top_v_eq) > crite_eq:
#                                            print('EQ:')
#                                        tr.plot(type='relative',color='g')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                                        event_day.append(tr.stats.starttime)
                                        
                                        on_off_all[event_count][0]=on_off1[x,0]
                                        on_off_all[event_count][1]=on_off1[x,1]
                                        event_count += 1
                                        on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))


#                                        file_sta.write(str(tr.stats.starttime + 8))
#                                        file_sta.write("\n")
#                                        event_eq.append(tr.stats.starttime) 
#                                        count+=1

    
#print(event)    
#file_sta.close()
#print(count)

    #%% LB02 only events
     
    LB02_start=[]
    LB02_end=[]                         
    for q in range(0,len(on_off2)):
        
        value = on_off2[q,0]       
        near,ind=find_nearest(on_off_all[:,0], value)

        if abs(value-near) > 100*15: #15s          

            LB02_start.append(on_off2[q,0]/sr)
            LB02_end.append(on_off2[q,1]/sr)

    for x in range(0,len(LB02_start)):
        tr = trace2.slice(starttime=start+(LB02_start[x]) , endtime=start+(LB02_end[x]))
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
                              
                                top_v_eq,top_eq,corell_eq = corel(trc2_e,trs_e,shift)
                                
                                top_v_n,top_n,corell_n = corel(trc3_e,trs_e,shift)
                           
                                if abs(top_v) < crite:
                                    if abs(top_v_n) > crite_n:
                                        if abs(top_v_eq) > crite_eq:
    #                                            print('EQ:')
    #                                        print('cf',cf,'peak',peak,'bwid',bwid,'Xc',top_v_eq,'rat',amp_rat)
    #                                        tr.plot(type='relative',color='g')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                                            event_day.append(tr.stats.starttime)
                                            
                                            on_off_all[event_count][0]=LB02_start[x]
                                            on_off_all[event_count][1]=LB02_end[x]
                                            event_count += 1
                                            on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))

#                                            print('lb03 event added')
        except:
            print('find at edge of data')               

    #%% LB03 only events
     
    LB03_start=[]
    LB03_end=[]                         
    for q in range(0,len(on_off3)):
        
        value = on_off3[q,0]
        near,ind=find_nearest(on_off_all[:,0], value)

        if abs(value-near) > 100*15: #15s          

            LB03_start.append(on_off3[q,0]/sr)
            LB03_end.append(on_off3[q,1]/sr)

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
                              
                                top_v_eq,top_eq,corell_eq = corel(trc2_e,trs_e,shift)
                                
                                top_v_n,top_n,corell_n = corel(trc3_e,trs_e,shift)
                           
                                if abs(top_v) < crite:
                                    if abs(top_v_n) > crite_n:
                                        if abs(top_v_eq) > crite_eq:
    #                                            print('EQ:')
    #                                        print('cf',cf,'peak',peak,'bwid',bwid,'Xc',top_v_eq,'rat',amp_rat)
    #                                        tr.plot(type='relative',color='g')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                                            event_day.append(tr.stats.starttime)
                                            
                                            on_off_all[event_count][0]=LB03_start[x]
                                            on_off_all[event_count][1]=LB03_end[x]
                                            event_count += 1
                                            on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))

#                                            print('lb03 event added')
        except:
            print('find at edge of data')
            
    #%% LB04 only events
     
    LB04_start=[]
    LB04_end=[]                         
    for q in range(0,len(on_off4)):
        
        value = on_off4[q,0]
        near,ind=find_nearest(on_off_all[:,0], value)

        if abs(value-near) > 100*15: #15s          

            LB04_start.append(on_off4[q,0]/sr)
            LB04_end.append(on_off4[q,1]/sr)

    for x in range(0,len(LB04_start)):
        tr = trace4.slice(starttime=start+(LB04_start[x]) , endtime=start+(LB04_end[x]))
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
                              
                                top_v_eq,top_eq,corell_eq = corel(trc2_e,trs_e,shift)
                                
                                top_v_n,top_n,corell_n = corel(trc3_e,trs_e,shift)
                           
                                if abs(top_v) < crite:
                                    if abs(top_v_n) > crite_n:
                                        if abs(top_v_eq) > crite_eq:
    #                                            print('EQ:')
    #                                        print('cf',cf,'peak',peak,'bwid',bwid,'Xc',top_v_eq,'rat',amp_rat)
    #                                        tr.plot(type='relative',color='g')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                                            event_day.append(tr.stats.starttime)
                                            
                                            on_off_all[event_count][0]=LB04_start[x]
                                            on_off_all[event_count][1]=LB04_end[x]
                                            event_count += 1
                                            on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))

#                                            print('lb03 event added')
        except:
            print('find at edge of data')
            
    
    #%% LB05 only events
     
    LB05_start=[]
    LB05_end=[]                         
    for q in range(0,len(on_off5)):
        
        value = on_off5[q,0]
        near,ind=find_nearest(on_off_all[:,0], value)

        if abs(value-near) > 100*15: #15s          

            LB05_start.append(on_off5[q,0]/sr)
            LB05_end.append(on_off5[q,1]/sr)

    for x in range(0,len(LB05_start)):
        tr = trace5.slice(starttime=start+(LB05_start[x]) , endtime=start+(LB05_end[x]))
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
                              
                                top_v_eq,top_eq,corell_eq = corel(trc2_e,trs_e,shift)
                                
                                top_v_n,top_n,corell_n = corel(trc3_e,trs_e,shift)
                           
                                if abs(top_v) < crite:
                                    if abs(top_v_n) > crite_n:
                                        if abs(top_v_eq) > crite_eq:
    #                                            print('EQ:')
    #                                        print('cf',cf,'peak',peak,'bwid',bwid,'Xc',top_v_eq,'rat',amp_rat)
    #                                        tr.plot(type='relative',color='g')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                                            event_day.append(tr.stats.starttime)
                                            
                                            on_off_all[event_count][0]=LB05_start[x]
                                            on_off_all[event_count][1]=LB05_end[x]
                                            event_count += 1
                                            on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))

#                                            print('lb03 event added')
        except:
            print('find at edge of data')
            
    #%% LB06 only events
     
    LB06_start=[]
    LB06_end=[]                         
    for q in range(0,len(on_off6)):
        
        value = on_off6[q,0]
        near,ind=find_nearest(on_off_all[:,0], value)

        if abs(value-near) > 100*15: #15s          

            LB06_start.append(on_off6[q,0]/sr)
            LB06_end.append(on_off6[q,1]/sr)

    for x in range(0,len(LB06_start)):
        tr = trace6.slice(starttime=start+(LB06_start[x]) , endtime=start+(LB06_end[x]))
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
                              
                                top_v_eq,top_eq,corell_eq = corel(trc2_e,trs_e,shift)
                                
                                top_v_n,top_n,corell_n = corel(trc3_e,trs_e,shift)
                           
                                if abs(top_v) < crite:
                                    if abs(top_v_n) > crite_n:
                                        if abs(top_v_eq) > crite_eq:
    #                                            print('EQ:')
    #                                        print('cf',cf,'peak',peak,'bwid',bwid,'Xc',top_v_eq,'rat',amp_rat)
    #                                        tr.plot(type='relative',color='g')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                                            event_day.append(tr.stats.starttime)
                                            
                                            on_off_all[event_count][0]=LB06_start[x]
                                            on_off_all[event_count][1]=LB06_end[x]
                                            event_count += 1
                                            on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))

#                                            print('lb03 event added')
        except:
            print('find at edge of data')
            
    
    #%% LB07 only events
     
    LB07_start=[]
    LB07_end=[]                         
    for q in range(0,len(on_off7)):
        
        value = on_off7[q,0]
        near,ind=find_nearest(on_off_all[:,0], value)

        if abs(value-near) > 100*15: #15s          

            LB07_start.append(on_off7[q,0]/sr)
            LB07_end.append(on_off7[q,1]/sr)

    for x in range(0,len(LB07_start)):
        tr = trace7.slice(starttime=start+(LB07_start[x]) , endtime=start+(LB07_end[x]))
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
                              
                                top_v_eq,top_eq,corell_eq = corel(trc2_e,trs_e,shift)
                                
                                top_v_n,top_n,corell_n = corel(trc3_e,trs_e,shift)
                           
                                if abs(top_v) < crite:
                                    if abs(top_v_n) > crite_n:
                                        if abs(top_v_eq) > crite_eq:
    #                                            print('EQ:')
    #                                        print('cf',cf,'peak',peak,'bwid',bwid,'Xc',top_v_eq,'rat',amp_rat)
    #                                        tr.plot(type='relative',color='g')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                                            event_day.append(tr.stats.starttime)
                                            
                                            on_off_all[event_count][0]=LB07_start[x]
                                            on_off_all[event_count][1]=LB07_end[x]
                                            event_count += 1
                                            on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))

#                                            print('lb03 event added')
        except:
            print('find at edge of data')
            
    #%% LB08 only events
     
    LB08_start=[]
    LB08_end=[]                         
    for q in range(0,len(on_off8)):
        
        value = on_off8[q,0]
        near,ind=find_nearest(on_off_all[:,0], value)

        if abs(value-near) > 100*15: #15s          

            LB08_start.append(on_off8[q,0]/sr)
            LB08_end.append(on_off8[q,1]/sr)

    for x in range(0,len(LB08_start)):
        tr = trace8.slice(starttime=start+(LB08_start[x]) , endtime=start+(LB08_end[x]))
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
                              
                                top_v_eq,top_eq,corell_eq = corel(trc2_e,trs_e,shift)
                                
                                top_v_n,top_n,corell_n = corel(trc3_e,trs_e,shift)
                           
                                if abs(top_v) < crite:
                                    if abs(top_v_n) > crite_n:
                                        if abs(top_v_eq) > crite_eq:
    #                                            print('EQ:')
    #                                        print('cf',cf,'peak',peak,'bwid',bwid,'Xc',top_v_eq,'rat',amp_rat)
    #                                        tr.plot(type='relative',color='g')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                                            event_day.append(tr.stats.starttime)
                                            
                                            on_off_all[event_count][0]=LB08_start[x]
                                            on_off_all[event_count][1]=LB08_end[x]
                                            event_count += 1
                                            on_off_all = np.lib.pad(on_off_all, ((0,1),(0,0)), 'constant', constant_values=(0))

#                                            print('lb03 event added')
        except:
            print('find at edge of data')
            

  
            
#%%            
    event_day.sort()

    for p in range(0,len(event_day)):
        file_sta.write(str(event_day[p]))
        file_sta.write("\n")
file_sta.close()
print('end of scan')
        

        

    
