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
bwid_max=0
bwid_min=10000
band_start=0
band_end=0
bwid_tot=0
shift=200
step=2
crite_eq=0.19
crite=0.55
crite_n=0.2
#s_window=7*60*60 
window=30
fmin=3
fmax=12
#file_sta = open('/Users/william/Documents/scanner/sta_lta_RF_test.txt','w')
#%% Read in waveforms   
#RF 
sample = read('/Users/william/Documents/lb01/14_336z.mseed')
trace=sample[0]
trace.detrend(type='demean')
trace.filter("bandpass", freqmin=fmin,freqmax=fmax)

startw= sample[0].stats.starttime + 14*60*60 + 11*60 + 34 #time window start - several s before Pwave arrival
endw=startw + window
trc = trace.slice(starttime = startw  , endtime= endw) #cut out sample waveform with same window length as chosen event
trc_e = obspy.signal.filter.envelope(trc.data)
print('reference RF waveform')
trc.plot(type='relative',color='b', starttime=startw , endtime=endw)

#EQ
sample2 = read('/Users/william/Documents/lb01/14_328z.mseed')
trace2=sample2[0]
trace2.detrend(type='demean')
trace2.filter("bandpass", freqmin=fmin,freqmax=fmax)

start2= sample2[0].stats.starttime + 8*60*60 + 30*60 +10 #time window start - 2s before Pwave arrival
end2=start2 + window +10
trs2 = trace2.slice(starttime = start2  , endtime= end2) #cut out sample waveform with same window length as chosen event
trs2_e = obspy.signal.filter.envelope(trs2.data)
print('reference EQ waveform')
trs2.plot(type='relative',color='b', starttime=start2 , endtime=end2)

#noise
sample3 = read('/Users/william/Documents/lb01/14_335z.mseed')
trace3=sample3[0]
trace3.detrend(type='demean')
trace3.filter("bandpass", freqmin=fmin,freqmax=fmax)

start3= sample3[0].stats.starttime + 20*60*60 + 57*60 + 50 #time window start - several s before Pwave arrival
end3=start3 + window
trs3 = trace3.slice(starttime = start3  , endtime= end3) #cut out sample waveform with same window length as chosen event
trs3_e = obspy.signal.filter.envelope(trs3.data)
print('reference noise waveform')
trs3.plot(type='relative',color='b', starttime=start3 , endtime=end3)
#%%  
stream = sample.copy()
stream = stream.clear()

for file in glob.glob('/Users/william/Documents/lb01/14_33*.mseed'):

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
trig_on=12                                            #8
trig_off=0.2                                         #0.2

event=[]
event_eq=[]
for t in range(0,len(stream),1): 
    print('Day',t+1,'of',len(stream))
#    file_sta.write("\n")
    trs=stream[t]
    trace=trs
    
    #window endpoints
    start= stream[t].stats.starttime   #time window start 
#    end=stream[t].stats.endtime 
#    trs = trace.slice(starttime = start  , endtime= end) #cut out sample waveform with same window length as chosen event
#    trs_e = obspy.signal.filter.envelope(trs.data)
    #print('reference waveform')
    
    
#    trs.plot(type='relative',color='b')#, starttime=start , endtime=end)
    st=trs.data
    cft=recursive_sta_lta(st, nsta, nlta)
                                       
#    plot_trigger(trs, cft, trig_on, trig_off) 
    
    on_off = trigger_onset(cft,trig_on,trig_off)
    
    for x in range(0,len(on_off)):
        tr = trace.slice(starttime=start+(on_off[x,0]/sr)-8 , endtime=start+(on_off[x,1]/sr)+8)
        trs_e = obspy.signal.filter.envelope(tr.data)
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
        
        # bandwidth 50%
        values = famp
        searchval = max(famp)/2
        ii = np.where(values >= searchval)[0]
        
        interval=tr_len
        above=[]
        for xx in range(0,len(ii)):
            if ii[xx]/interval < 50:
                above.append(ii[xx]/interval)
        bwid=max(above)-min(above)
        
        band_start+= min(above)
        band_end+= max(above)
        if bwid > bwid_max:
            bwid_max=bwid
        if bwid < bwid_min:
            bwid_min=bwid
        bwid_tot+= bwid
                
        
        tr.detrend(type='demean')
        amp=abs(tr.max())
        whole = sum(abs(trs_e))
        av_amp = whole/(100*tr_len)
        amp_rat = amp/av_amp
        # only look for things above noise level 
        
        if amp > 600:
            if amp_rat < 15:
                if 5 < cf < 10 :  
                    if 2 < peak < 14:
                        if 1 < bwid < 13:
                            corel=correlate(trc_e, trs_e, shift, demean=True,normalize='naive',domain='time' )
                            top=corel.argmax()
                            bot=corel.argmin()
                            top_v = corel[top]
                            bot_v = corel[bot]  
                            
                            if abs(bot_v) > abs(top_v):
                                top_v=abs(bot_v)
                                top=100-bot #goes from -100 to +100 but list goes from 0 to 200
                                #print('data-eq:', top_v)
                            
                            
                            corel_eq=correlate(trs2_e, trs_e, shift, demean=True,normalize='naive',domain='time' )
                            top_eq=corel_eq.argmax()
                            bot_eq=corel_eq.argmin()
                            top_v_eq = corel_eq[top_eq]
                            bot_v_eq = corel_eq[bot_eq]  
                            
                            if abs(bot_v_eq) > abs(top_v_eq):
                                top_v_eq=abs(bot_v_eq)
                                top_eq=100-bot_eq #goes from -100 to +100 but list goes from 0 to 200
                                #print('data-eq:', top_v)
                                
                            corel_n=correlate(trs3_e, trs_e, shift, demean=True,normalize='naive',domain='time' )
                            top_n=corel_n.argmax()
                            bot_n=corel_n.argmin()
                            top_v_n = corel_eq[top_n]
                            bot_v_n = corel_eq[bot_n]  
                            
                            if abs(bot_v_n) > abs(top_v_n):
                                top_v_n=abs(bot_v_n)
                                top_n=100-bot_n #goes from -100 to +100 but list goes from 0 to 200
                                #print('data-eq:', top_v)
                                
                                
#                            print(cf,peak,bwid,top_v,top_v_eq,top_v_n,amp_rat)
                            
#                             
#                            if 5 < cf < 8.2:
#                                if 1 < bwid < 6:
#                                    if 3.2 < peak < 8.5:
#                                        if abs(top_v_n) > crite_n:
#                                            if abs(top_v) > crite:
#                                                print('RF:')
#                    
#                                                tr.plot(type='relative',color='b')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
#                                                event.append(tr.stats.starttime)
    #                                    else:
    #                                        print('noise: low crite')
    #                                        tr.plot(type='relative',color='r')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
    #                                else:
    #                                    print('noise: wrong peak')
    #                                    tr.plot(type='relative',color='r')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
    #                            else:
    #                                print('noise: wrong bwid')
    #                                tr.plot(type='relative',color='r')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
    #                            
    #                        else:
    #                            print('noise: wrong cf')
    #                            tr.plot(type='relative',color='r')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
    
    
                            
                            if 3 < bwid < 13:
                                if abs(top_v) < crite:
                                    if abs(top_v_n) > crite_n:
                                        if abs(top_v_eq) > crite_eq:
#                                            print('EQ:')
                                            print('cf',cf,'peak',peak,'bwid',bwid,'Xc',top_v_eq,'rat',amp_rat)
                                            tr.plot(type='relative',color='g')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
#                                            event.append(tr.stats.starttime)
    #                            else:
    #                                print('noise: low crite')
    #                                tr.plot(type='relative',color='r')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
    #                        else:
    #                            print('noise: wrong bwid')
    #                            tr.plot(type='relative',color='r')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
                                            
                                            
                                            
                            if abs(top_v) < crite:
                                if abs(top_v_eq) < crite_eq: 
                                    print('noise:')
                
#                                    tr.plot(type='relative',color='r')#, starttime=start+(on_off[x,0]/sr)-10 , endtime=start+(on_off[x,1]/sr)+10)
#                                    event.append(tr.stats.starttime)
    
#print(event)    
#file_sta.close()
    
    
    # DO NOT DELETE
