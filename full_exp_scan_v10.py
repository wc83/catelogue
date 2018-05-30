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
from obspy import UTCDateTime

T1= time.clock()

#%%
file_exp = open('/Users/william/Documents/scanner/all_stations/full_exp_scan_v10.txt','w')
file_cor2 = open('/Users/william/Documents/scanner/all_stations/full_correlated2_exp_scan_v10.txt','w')
file_lb01= open('/Users/william/Documents/scanner/LB01/LB01_full_exp_scan_v10.txt','w')
file_lb02= open('/Users/william/Documents/scanner/LB02/LB02_full_exp_scan_v10.txt','w')
file_lb03= open('/Users/william/Documents/scanner/LB03/LB03_full_exp_scan_v10.txt','w')
file_lb04= open('/Users/william/Documents/scanner/LB04/LB04_full_exp_scan_v10.txt','w')
file_lb05= open('/Users/william/Documents/scanner/LB05/LB05_full_exp_scan_v10.txt','w')
file_lb06= open('/Users/william/Documents/scanner/LB06/LB06_full_exp_scan_v10.txt','w')


   #%% constants    
shift=100
step=2
crite=0.5 
window=50
fmin=0.1
fmax=10
#%% Reference stacked waveforms     

#trs_e1,trs_e2,trs_e3,trs_e4,trs_e5,trs_e6,trs_e7=get_reference(fmin,fmax)
trs_e1,trs_e2,trs_e3,trs_e4,trs_e5,trs_e6=stacked_LBz_exp() 


#%% open empty streams for each station
sample = read('/Users/william/Documents/lb01/14_336z.mseed')
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

#%% Read in Each station
num_active=[]
for x in range(14,15):
    
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

#%% Open matracies

Events_per_day=[]
Event_store=[]
All_events=[]

#%% Determine all stations as inactive until data is detected
switch1=0
switch2=0
switch3=0
switch4=0
switch5=0
switch6=0

#%% Loop over whole time one day at a time
cor_events=np.zeros(shape=(1,14))
cor_count=-1
data_array=[]
for p in range(0,len(stream1)):
    file_exp.write("\n")
    file_cor2.write("\n")
    print('Day',p+1,'of',len(stream1))
    Events_today=np.zeros(shape=(1,8))
    todays_event=[]  
    correlated_today=[]                                    
    
    Number_of_events_today=0
    find_count=0    
    
    trace1=stream1[p]
    trace2=stream2[p]
    trace3=stream3[p]
    trace4=stream4[p]
    trace5=stream5[p]
    trace6=stream6[p]

#%% Station LB01 scan    
    if trace1.stats.npts > 201: 
        trace1.detrend(type='linear')
        trace1.detrend(type='demean')
        trace1.filter("bandpass", freqmin=fmin,freqmax=fmax)
        if switch1 == 0:
            sw_start1=trace1.stats.starttime
            tr_av1 = trace1.slice(starttime = sw_start1  , endtime= sw_start1 +30*60)
            tr_av1.detrend(type='linear')
            tr_av1.detrend(type='demean')
            av_amp1=abs(sum(abs(tr_av1.data))/len(tr_av1))
            av_amp1=av_amp1*5
            if av_amp1 > 100:
                switch1=1
        if switch1 ==1:
        
            begin = trace1.stats.starttime
            add=0
            lastcorr=0
            expt=0
            dl=trace1.stats.endtime - trace1.stats.starttime 
            # Loop over single day, scanning for waveforms
            for t in range(0,int(trace1.stats.npts/100), step):
                view_start = begin + t + add
                view_end = view_start + window
                if view_start > trace1.stats.endtime - window :
                    break
        
                trc = trace1.slice(starttime = view_start  , endtime= view_end)
                trc_save = trace1.slice(starttime = view_start -10 , endtime= view_end +10)
                trc.detrend(type='linear')
                trc.detrend(type='demean')
                trc_save.detrend(type='linear')
                trc_save.detrend(type='demean')
                # only look for things above noise level 
                amp=abs(trc.max())
                
                if amp > av_amp1:                
                    peak,cf,bwid,bwid25 = freq_info(trc.data,trc.stats.starttime,trc.stats.endtime)         
                    if (0.75 < cf < 2.75) and (0.2 < peak < 2.5) and (0.3 < bwid < 2.5): 
                             #correlate between data and earthquake envelopes 
                        trc_e = obspy.signal.filter.envelope(trc.data) 
                        top_v,top,corell = corel(trs_e1,trc_e,shift)
                        
                        if 0.0 in corell : # in the case of missing data - which will cause function to crash
                            break   
                        
                        expt_new = view_start - (top/100) +10

                        if abs(top_v) > crite and top_v > lastcorr:
                                lastcorr = top_v
                                expt = expt_new                                    
                            
                        if lastcorr > top_v:

                            rt=expt.timestamp	
                            jd=expt.julday
                            yr=expt.year
                            mo=expt.month
                            da=expt.day
                            hr=expt.hour
                            mi=expt.minute
                            se=expt.second
                            ms=expt.microsecond/10000
                            row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,x])
                            file_lb01.write(str(expt))
                            file_lb01.write("\n")
                            
                            
                            todays_event.append(expt)
                            Event_store.append(row)
                            Events_today[find_count][0]=rt
                            Events_today[find_count][1]=p
                            Events_today[find_count][2]=1
                            Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                            All_events.append(expt)
                            find_count += 1
                            Number_of_events_today += 1
                            
                            if num_active[p] == 1:
                                cor_events[cor_count][0]=rt
                                cor_events[cor_count][1]="%1d" % p
                                cor_events[cor_count][2]="%1d" % yr
                                cor_events[cor_count][3]="%1d" % jd
                                cor_events[cor_count][4]="%1d" % hr
                                cor_events[cor_count][5]="%1d" % mi
                                cor_events[cor_count][6]="%1d" % se
                                cor_events[cor_count][7]="%1d" % ms
                                cor_events[cor_count][8]="%1d" % 1
                                cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
                                correlated_today.append(expt)
#                                        data_array.append(trc_save)
                                cor_count += 1
                            lastcorr=0
                            add += 60-step               #skip 60s to avoid double catch

                else:
                    add += 38
                    lastcorr=0 
    else:
        switch1=0
 #%% Station 2 scan  
    if trace2.stats.npts > 201:
        trace2.detrend(type='linear')
        trace2.detrend(type='demean')
        trace2.filter("bandpass", freqmin=fmin,freqmax=fmax)
        if switch2 == 0:
            sw_start2=trace2.stats.starttime
            tr_av2 = trace2.slice(starttime = sw_start2  , endtime= sw_start2 +30*60)
            tr_av2.detrend(type='linear')
            tr_av2.detrend(type='demean')
            av_amp2=abs(sum(abs(tr_av2.data))/len(tr_av2))
            av_amp2=av_amp2*5
            if av_amp2 > 100:
                switch2=1
        if switch2 ==1:
            
            begin = trace2.stats.starttime
            add=0
            lastcorr=0     
            expt=0
            dl=trace2.stats.endtime - trace2.stats.starttime     
            for t in range(0,int(trace2.stats.npts/100), step):
                view_start = begin + t + add
                view_end = view_start + window
                if view_start > trace2.stats.endtime - window :
                    break
        
                trc2 = trace2.slice(starttime = view_start  , endtime= view_end)
                trc2_save = trace1.slice(starttime = view_start -10 , endtime= view_end +10)
                trc2.detrend(type='linear')
                trc2.detrend(type='demean')
                trc2_save.detrend(type='linear')
                trc2_save.detrend(type='demean')
                # only look for things above noise level 
                amp=abs(trc2.max())
                
                if amp > av_amp2:  
                    peak2,cf2,bwid2,bwid252 = freq_info(trc2.data,trc2.stats.starttime,trc2.stats.endtime)      
                    if (0.75 < cf2 < 2.75) and (0.4 < peak2 < 3) and (0.2 < bwid2 < 2.5):
                             #correlate between data and earthquake envelopes 
                        trc2_e = obspy.signal.filter.envelope(trc2.data) 
                        top_v,top,corell = corel(trs_e2,trc2_e,shift)
                        
                        if 0.0 in corell : # in the case of missing data - which will cause function to crash
                            break            
                        
                        expt_new = view_start - (top/100) +10
            #            if the correlation is positive, exp found
                        if abs(top_v) > crite and top_v > lastcorr:
                            lastcorr = top_v
                            expt = expt_new  

                        if lastcorr > top_v:
                        
                            file_lb02.write(str(expt))
                            file_lb02.write("\n")
#                                
                            rt=expt.timestamp	
                            jd=expt.julday
                            yr=expt.year
                            mo=expt.month
                            da=expt.day
                            hr=expt.hour
                            mi=expt.minute
                            se=expt.second
                            ms=expt.microsecond/10000	
                            near,ind=find_nearest(Events_today[:,0], rt)
                            if abs(rt-near) > 60: #60s
  
                                row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,x])
                                
                                
                                todays_event.append(expt)
                                Event_store.append(row)
                                Events_today[find_count][0]=rt
                                Events_today[find_count][1]=p
                                Events_today[find_count][3]=1
                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                All_events.append(expt)
                                find_count += 1
                                Number_of_events_today += 1
                            
                            else:
                                Events_today[ind][3]=1
                            
                            if num_active[p] > 1:
                                if abs(rt-near) <= 6:
                                    near_cor,inde=find_nearest(cor_events[:,0], rt)
                                    if abs(rt-near_cor) >  30:
                                        cor_events[cor_count][0]=rt
                                        cor_events[cor_count][1]="%1d" % p
                                        cor_events[cor_count][2]="%1d" % yr
                                        cor_events[cor_count][3]="%1d" % jd
                                        cor_events[cor_count][4]="%1d" % hr
                                        cor_events[cor_count][5]="%1d" % mi
                                        cor_events[cor_count][6]="%1d" % se
                                        cor_events[cor_count][7]="%1d" % ms
                                        cor_events[cor_count][8]="%1d" % Events_today[ind][2]
                                        cor_events[cor_count][9]="%1d" % Events_today[ind][3]
                                        cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        correlated_today.append(expt)
#                                                data_array.append(trc2_save)
                                        cor_count += 1
                                        
                            else:
                                cor_events[cor_count][0]=rt
                                cor_events[cor_count][1]="%1d" % p
                                cor_events[cor_count][2]="%1d" % yr
                                cor_events[cor_count][3]="%1d" % jd
                                cor_events[cor_count][4]="%1d" % hr
                                cor_events[cor_count][5]="%1d" % mi
                                cor_events[cor_count][6]="%1d" % se
                                cor_events[cor_count][7]="%1d" % ms
                                cor_events[cor_count][9]="%1d" % 1
                                cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
                                correlated_today.append(expt)
#                                        data_array.append(trc2_save)
                                cor_count += 1
                            
                            lastcorr =0
                            add += 60-step               #skip 60s to avoid double catch

                else:
                    add += 38
                    lastcorr=0 
    else:
        switch2=0
 
 #%% Station 3 scan  
    if trace3.stats.npts > 201:
        trace3.detrend(type='linear')
        trace3.detrend(type='demean')
        trace3.filter("bandpass", freqmin=fmin,freqmax=fmax)
        if switch3 == 0:
            sw_start3=trace3.stats.starttime
            tr_av3 = trace3.slice(starttime = sw_start3  , endtime= sw_start3 +30*60)
            tr_av3.detrend(type='linear')
            tr_av3.detrend(type='demean')
            av_amp3=abs(sum(abs(tr_av3.data))/len(tr_av3))
            av_amp3=av_amp3*5
            if av_amp3 > 100:
                switch3=1
        if switch3 ==1:
        
            begin = trace3.stats.starttime
            add=0
            lastcorr=0  
            expt=0
            dl=trace3.stats.endtime - trace3.stats.starttime     
            for t in range(0,int(trace3.stats.npts/100), step):
                view_start = begin + t + add
                view_end = view_start + window
                if view_start > trace3.stats.endtime - window :
                    break
        
                trc3 = trace3.slice(starttime = view_start  , endtime= view_end)
                trc3_save = trace1.slice(starttime = view_start -10 , endtime= view_end +10)
                trc3.detrend(type='linear')
                trc3.detrend(type='demean')
                trc3_save.detrend(type='linear')
                trc3_save.detrend(type='demean')
                # only look for things above noise level 
                amp=abs(trc3.max())
                
                if amp > av_amp3:       
                    peak3,cf3,bwid3,bwid253 = freq_info(trc3.data,trc3.stats.starttime,trc3.stats.endtime)      
                    if (0.9 < cf3 < 4) and (0.2 < peak3 < 2.5) and (0.2 < bwid3 < 2.75):
                             #correlate between data and earthquake envelopes 
                        trc3_e = obspy.signal.filter.envelope(trc3.data) 
                        top_v,top,corell = corel(trs_e3,trc3_e,shift)
                        
                        if 0.0 in corell : # in the case of missing data - which will cause function to crash
                            break            
 
                        expt_new = view_start - (top/100) +10

                        if abs(top_v) > crite and top_v > lastcorr:
                                lastcorr = top_v
                                expt = expt_new                                    
                            
                        if lastcorr > top_v:
                            
                            file_lb03.write(str(expt))
                            file_lb03.write("\n")
                            
                            rt=expt.timestamp	
                            jd=expt.julday
                            yr=expt.year
                            mo=expt.month
                            da=expt.day
                            hr=expt.hour
                            mi=expt.minute
                            se=expt.second
                            ms=expt.microsecond/10000	
                            near,ind=find_nearest(Events_today[:,0], rt)
                            if abs(rt-near) > 60: #60s
  
                                row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,x])
                                
                                
                                todays_event.append(expt)
                                Event_store.append(row)
                                Events_today[find_count][0]=rt
                                Events_today[find_count][1]=p
                                Events_today[find_count][4]=1
                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                All_events.append(expt)
                                find_count += 1
                                Number_of_events_today += 1
                            
                            else:
                                Events_today[ind][4]=1

                            if num_active[p] > 1:
                                if abs(rt-near) <= 6:
                                    near_cor,inde=find_nearest(cor_events[:,0], rt)
                                    if abs(rt-near_cor) >  30:
                                        cor_events[cor_count][0]=rt
                                        cor_events[cor_count][1]="%1d" % p
                                        cor_events[cor_count][2]="%1d" % yr
                                        cor_events[cor_count][3]="%1d" % jd
                                        cor_events[cor_count][4]="%1d" % hr
                                        cor_events[cor_count][5]="%1d" % mi
                                        cor_events[cor_count][6]="%1d" % se
                                        cor_events[cor_count][7]="%1d" % ms
                                        cor_events[cor_count][8]="%1d" % Events_today[ind][2]
                                        cor_events[cor_count][9]="%1d" % Events_today[ind][3]
                                        cor_events[cor_count][10]="%1d" % Events_today[ind][4]
                                        cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        correlated_today.append(expt)
#                                                data_array.append(trc3_save)
                                        cor_count += 1
                                    else:
                                        cor_events[inde][10]="%1d" % 1
#                                                data_array.append(trc3_save)
                                        
                                       
                            else:
                                cor_events[cor_count][0]=rt
                                cor_events[cor_count][1]="%1d" % p
                                cor_events[cor_count][2]="%1d" % yr
                                cor_events[cor_count][3]="%1d" % jd
                                cor_events[cor_count][4]="%1d" % hr
                                cor_events[cor_count][5]="%1d" % mi
                                cor_events[cor_count][6]="%1d" % se
                                cor_events[cor_count][7]="%1d" % ms
                                cor_events[cor_count][10]="%1d" % 1
                                cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
                                correlated_today.append(expt)
#                                        data_array.append(trc3_save)
                                cor_count += 1
                            lastcorr=0    
                            add += 60-step               #skip 60s to avoid double catch

                else:
                    add += 38
                    lastcorr=0 
    else:
        switch3=0
                
                 #%% Station 4 scan  
    if trace4.stats.npts > 201:
        trace4.detrend(type='linear')
        trace4.detrend(type='demean')
        trace4.filter("bandpass", freqmin=fmin,freqmax=fmax)
        if switch4 == 0:
            sw_start4=trace4.stats.starttime
            tr_av4 = trace4.slice(starttime = sw_start4  , endtime= sw_start4 +30*60)
            tr_av4.detrend(type='linear')
            tr_av4.detrend(type='demean')
            av_amp4=abs(sum(abs(tr_av4.data))/len(tr_av4))
            av_amp4=av_amp4*5
            if av_amp4 > 100:
                switch4=1
        if switch4 ==1:
        
            begin = trace4.stats.starttime
            add=0
            lastcorr=0  
            expt=0
            dl=trace4.stats.endtime - trace4.stats.starttime     
            for t in range(0,int(trace4.stats.npts/100), step):
                view_start = begin + t + add
                view_end = view_start + window
                if view_start > trace4.stats.endtime - window :
                    break
        
                trc4 = trace4.slice(starttime = view_start  , endtime= view_end)
                trc4_save = trace1.slice(starttime = view_start -10 , endtime= view_end +10)
                trc4.detrend(type='linear')
                trc4.detrend(type='demean')
                trc4_save.detrend(type='linear')
                trc4_save.detrend(type='demean')
                # only look for things above noise level 
                amp=abs(trc4.max())
                
                if amp > av_amp4:       
                    peak4,cf4,bwid4,bwid254 = freq_info(trc4.data,trc4.stats.starttime,trc4.stats.endtime)      
                    if (0.9 < cf4 < 5) and (0.4 < peak4 < 3.5) and (0.3 < bwid4 < 4): 
                             #correlate between data and earthquake envelopes 
                        trc4_e = obspy.signal.filter.envelope(trc4.data) 
                        top_v,top,corell = corel(trs_e4,trc4_e,shift)
                       
                        if 0.0 in corell : # in the case of missing data - which will cause function to crash
                            break            
 
                        expt_new = view_start - (top/100) +10

                        if abs(top_v) > crite and top_v > lastcorr:
                                lastcorr = top_v
                                expt = expt_new                                    
                            
                        if lastcorr > top_v:

                            file_lb04.write(str(expt))
                            file_lb04.write("\n")
                            
                            rt=expt.timestamp	
                            jd=expt.julday
                            yr=expt.year
                            mo=expt.month
                            da=expt.day
                            hr=expt.hour
                            mi=expt.minute
                            se=expt.second
                            ms=expt.microsecond/10000
                            
                            near,ind=find_nearest(Events_today[:,0], rt)
                            if abs(rt-near) > 60: #60s
  
                                row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,x])
                                
                                todays_event.append(expt)
                                Event_store.append(row)
                                Events_today[find_count][0]=rt
                                Events_today[find_count][1]=p
                                Events_today[find_count][5]=1
                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                All_events.append(expt)
                                find_count += 1
                                Number_of_events_today += 1
                            else:
                                Events_today[ind][5]=1
                            
                            if num_active[p] > 1:
                                if abs(rt-near) <= 6:
                                    near_cor,inde=find_nearest(cor_events[:,0], rt)
                                    if abs(rt-near_cor) >  30:
                                        cor_events[cor_count][0]=rt
                                        cor_events[cor_count][1]="%1d" % p
                                        cor_events[cor_count][2]="%1d" % yr
                                        cor_events[cor_count][3]="%1d" % jd
                                        cor_events[cor_count][4]="%1d" % hr
                                        cor_events[cor_count][5]="%1d" % mi
                                        cor_events[cor_count][6]="%1d" % se
                                        cor_events[cor_count][7]="%1d" % ms
                                        cor_events[cor_count][8]="%1d" % Events_today[ind][2]
                                        cor_events[cor_count][9]="%1d" % Events_today[ind][3]
                                        cor_events[cor_count][10]="%1d" % Events_today[ind][4]
                                        cor_events[cor_count][11]="%1d" % Events_today[ind][5]
                                        cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        correlated_today.append(expt)
#                                                data_array.append(trc4_save)
                                        cor_count += 1
                                    else:
                                        cor_events[inde][11]="%1d" % 1
#                                                data_array.append(trc4_save)
                            else:
                                cor_events[cor_count][0]=rt
                                cor_events[cor_count][1]="%1d" % p
                                cor_events[cor_count][2]="%1d" % yr
                                cor_events[cor_count][3]="%1d" % jd
                                cor_events[cor_count][4]="%1d" % hr
                                cor_events[cor_count][5]="%1d" % mi
                                cor_events[cor_count][6]="%1d" % se
                                cor_events[cor_count][7]="%1d" % ms
                                cor_events[cor_count][11]="%1d" % 1
                                cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
                                correlated_today.append(expt)
#                                        data_array.append(trc4_save)
                                cor_count += 1
                                
                            lastcorr=0    
                            add += 60-step               #skip 60s to avoid double catch

                else:
                    add += 38
                    lastcorr=0 
    else:
        switch4=0
                
    #%% Station 5 scan  
    if trace5.stats.npts > 201:
        trace5.detrend(type='linear')
        trace5.detrend(type='demean')
        trace5.filter("bandpass", freqmin=fmin,freqmax=fmax)
        if switch5 == 0:
            sw_start5=trace5.stats.starttime
            tr_av5 = trace5.slice(starttime = sw_start5  , endtime= sw_start5 +30*60)
            tr_av5.detrend(type='linear')
            tr_av5.detrend(type='demean')
            av_amp5=abs(sum(abs(tr_av5.data))/len(tr_av5))
            av_amp5=av_amp5*5
            if av_amp5 > 100:
                switch5=1
        if switch5 ==1:
            
            begin = trace5.stats.starttime
            add=0
            lastcorr=0 
            expt=0
            dl=trace5.stats.endtime - trace5.stats.starttime    
            for t in range(0,int(trace5.stats.npts/100), step):
                view_start = begin + t + add
                view_end = view_start + window
                if view_start > trace5.stats.endtime - window :
                    break
        
                trc5 = trace5.slice(starttime = view_start  , endtime= view_end)
                trc5_save = trace1.slice(starttime = view_start -10 , endtime= view_end +10)
                trc5.detrend(type='linear')
                trc5.detrend(type='demean')
                trc5_save.detrend(type='linear')
                trc5_save.detrend(type='demean')
                # only look for things above noise level 
                amp=abs(trc5.max())
                
                if amp > av_amp5:       
                    peak5,cf5,bwid5,bwid255 = freq_info(trc5.data,trc5.stats.starttime,trc5.stats.endtime)      
                    if (0.9 < cf5 < 5) and (0.4 < peak5 < 3.5) and (0.3 < bwid5 < 3.5): 
                             #correlate between data and earthquake envelopes 
                        trc5_e = obspy.signal.filter.envelope(trc5.data) 
                        top_v,top,corell = corel(trs_e5,trc5_e,shift)
                        
                        if 0.0 in corell : # in the case of missing data - which will cause function to crash
                            break            
 
                        expt_new = view_start - (top/100) +10

                        if abs(top_v) > crite and top_v > lastcorr:
                                lastcorr = top_v
                                expt = expt_new                                    
                            
                        if lastcorr > top_v:

                            file_lb05.write(str(expt))
                            file_lb05.write("\n")
                            
                            rt=expt.timestamp	
                            jd=expt.julday
                            yr=expt.year
                            mo=expt.month
                            da=expt.day
                            hr=expt.hour
                            mi=expt.minute
                            se=expt.second
                            ms=expt.microsecond/10000
                            near,ind=find_nearest(Events_today[:,0], rt)
                            
                            if abs(rt-near) > 60: #60s
  
                                row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,x])
                                
                                todays_event.append(expt)
                                Event_store.append(row)
                                Events_today[find_count][0]=rt
                                Events_today[find_count][1]=p
                                Events_today[find_count][6]=1
                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                All_events.append(expt)
                                find_count += 1
                                Number_of_events_today += 1
                            else:
                                Events_today[ind][6]=1
                            
                            if num_active[p] > 1:
                                if abs(rt-near) <= 6:
                                    near_cor,inde=find_nearest(cor_events[:,0], rt)
                                    if abs(rt-near_cor) >  30:
                                        cor_events[cor_count][0]=rt
                                        cor_events[cor_count][1]="%1d" % p
                                        cor_events[cor_count][2]="%1d" % yr
                                        cor_events[cor_count][3]="%1d" % jd
                                        cor_events[cor_count][4]="%1d" % hr
                                        cor_events[cor_count][5]="%1d" % mi
                                        cor_events[cor_count][6]="%1d" % se
                                        cor_events[cor_count][7]="%1d" % ms
                                        cor_events[cor_count][8]="%1d" % Events_today[ind][2]
                                        cor_events[cor_count][9]="%1d" % Events_today[ind][3]
                                        cor_events[cor_count][10]="%1d" % Events_today[ind][4]
                                        cor_events[cor_count][11]="%1d" % Events_today[ind][5]
                                        cor_events[cor_count][12]="%1d" % Events_today[ind][6]
                                        cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        correlated_today.append(expt)
#                                                data_array.append(trc5_save)
                                        cor_count += 1
                                    else:
                                        cor_events[inde][12]="%1d" % 1
#                                                data_array.append(trc5_save)

                            else:
                                cor_events[cor_count][0]= rt
                                cor_events[cor_count][1]="%1d" % p
                                cor_events[cor_count][2]="%1d" % yr
                                cor_events[cor_count][3]="%1d" % jd
                                cor_events[cor_count][4]="%1d" % hr
                                cor_events[cor_count][5]="%1d" % mi
                                cor_events[cor_count][6]="%1d" % se
                                cor_events[cor_count][7]="%1d" % ms
                                cor_events[cor_count][12]="%1d" % 1
                                cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
                                correlated_today.append(expt)
#                                        data_array.append(trc5_save)
                                cor_count += 1
                            
                            lastcorr=0
                            add += 60-step               #skip 60s to avoid double catch

                else:
                    add += 38
                    lastcorr=0
    else:
        switch5=0
                
    #%% Station 6 scan  
    if trace6.stats.npts > 201:
        trace6.detrend(type='linear')
        trace6.detrend(type='demean')
        trace6.filter("bandpass", freqmin=fmin,freqmax=fmax)
        if switch6 == 0:
            sw_start6=trace6.stats.starttime
            tr_av6 = trace6.slice(starttime = sw_start6  , endtime= sw_start6 +30*60)
            tr_av6.detrend(type='linear')
            tr_av6.detrend(type='demean')
            av_amp6=abs(sum(abs(tr_av6.data))/len(tr_av6))
            av_amp6=av_amp6*5
            if av_amp6 > 100:
                switch6=1
        if switch6 ==1:
        
            begin = trace6.stats.starttime
            add=0
            lastcorr=0 
            expt=0
            dl=trace6.stats.endtime - trace6.stats.starttime    
            for t in range(0,int(trace6.stats.npts/100), step):
                view_start = begin + t + add
                view_end = view_start + window
                if view_start > trace6.stats.endtime - window :
                    break
        
                trc6 = trace6.slice(starttime = view_start  , endtime= view_end)
                trc6_save = trace1.slice(starttime = view_start -10 , endtime= view_end +10)
                trc6.detrend(type='linear')
                trc6.detrend(type='demean')
                trc6_save.detrend(type='linear')
                trc6_save.detrend(type='demean')
                # only look for things above noise level 
                amp=abs(trc6.max())
                
                if amp > av_amp6:       
                    peak6,cf6,bwid6,bwid256 = freq_info(trc6.data,trc6.stats.starttime,trc6.stats.endtime)      
                    if (0.65 < cf6 < 2.75) and (0.2 < peak6 < 2.5) and (0.3 < bwid6 < 3): 
                             #correlate between data and earthquake envelopes 
                        trc6_e = obspy.signal.filter.envelope(trc6.data) 
                        top_v,top,corell = corel(trs_e6,trc6_e,shift)
                        
                        if 0.0 in corell : # in the case of missing data - which will cause function to crash
                            break            
 
                        expt_new = view_start - (top/100) +10

                        if abs(top_v) > crite and top_v > lastcorr:
                                lastcorr = top_v
                                expt = expt_new                                    
                            
                        if lastcorr > top_v:

                            file_lb06.write(str(expt))
                            file_lb06.write("\n")
                            
                            rt=expt.timestamp	
                            jd=expt.julday
                            yr=expt.year
                            mo=expt.month
                            da=expt.day
                            hr=expt.hour
                            mi=expt.minute
                            se=expt.second
                            ms=expt.microsecond/10000
                            
                            near,ind=find_nearest(Events_today[:,0], rt)
                            if abs(rt-near) > 60: #60s
 
                                row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,x])
                                
                                todays_event.append(expt)
                                Event_store.append(row)
                                Events_today[find_count][0]=rt
                                Events_today[find_count][1]=p
                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                                All_events.append(expt)
                                find_count += 1
                                Number_of_events_today += 1
                            else:
                                Events_today[ind][7]=1
                            
                            if num_active[p] > 1:
                                if abs(rt-near) <= 6:
                                    near_cor,inde=find_nearest(cor_events[:,0], rt)
                                    if abs(rt-near_cor) >  30:
                                        cor_events[cor_count][0]=rt
                                        cor_events[cor_count][1]="%1d" % p
                                        cor_events[cor_count][2]="%1d" % yr
                                        cor_events[cor_count][3]="%1d" % jd
                                        cor_events[cor_count][4]="%1d" % hr
                                        cor_events[cor_count][5]="%1d" % mi
                                        cor_events[cor_count][6]="%1d" % se
                                        cor_events[cor_count][7]="%1d" % ms
                                        cor_events[cor_count][8]="%1d" % Events_today[ind][2]
                                        cor_events[cor_count][9]="%1d" % Events_today[ind][3]
                                        cor_events[cor_count][10]="%1d" % Events_today[ind][4]
                                        cor_events[cor_count][11]="%1d" % Events_today[ind][5]
                                        cor_events[cor_count][12]="%1d" % Events_today[ind][6]
                                        cor_events[cor_count][13]="%1d" % Events_today[ind][7]
                                        cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
                                        correlated_today.append(expt)
#                                                data_array.append(trc6_save)
                                        cor_count += 1
                                    else:
                                        cor_events[inde][13]="%1d" % 1
#                                                data_array.append(trc6_save)

                            else:
                                cor_events[cor_count][0]=rt
                                cor_events[cor_count][1]="%1d" % p
                                cor_events[cor_count][2]="%1d" % yr
                                cor_events[cor_count][3]="%1d" % jd
                                cor_events[cor_count][4]="%1d" % hr
                                cor_events[cor_count][5]="%1d" % mi
                                cor_events[cor_count][6]="%1d" % se
                                cor_events[cor_count][7]="%1d" % ms
                                cor_events[cor_count][13]="%1d" % 1
                                cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
                                correlated_today.append(expt)
#                                        data_array.append(trc6_save)
                                cor_count += 1
                            lastcorr=0                                        
                            add += 60-step               #skip 60s to avoid double catch
    
                else:
                    add += 38
                    lastcorr=0 
    else:
        switch6=0
                
#%% collect information from each station for single day and save in arrays                
    Events_per_day.append(Number_of_events_today)
    todays_event.sort()
    correlated_today.sort()
    for l in range(0,len(todays_event)):
        file_exp.write(str(todays_event[l]))
        file_exp.write("\n")
    
    if num_active[p] > 1:    
        for s in range(0,len(correlated_today)):
            file_cor2.write(str(correlated_today[s]))
            file_cor2.write("\n")       

#%%             
col=0       
cor_events=cor_events[np.argsort(cor_events[:,col])] 
if cor_events[0][0]== 0.0:
    cor_events = np.delete(cor_events, (0), axis=0)       
np.savetxt("/Users/william/Documents/scanner/all_stations/EXP_coincidence_V1.csv", cor_events,delimiter=",",header="Time_stamp,Day_number,Year,Day_of_year,Hour,Min,Sec,Milisec,LBO1,LB02,LB03,LB04,LB05,LB06")

for r in range(0,len(cor_events)):
    if cor_events[r][8]==1:
        tr1a=stream1[int(cor_events[r][1])]
        dstart=UTCDateTime(cor_events[r][0])
        tr1_save = tr1a.slice(starttime=dstart, endtime=dstart+80)
        tr1_save.detrend(type='demean')
        tr1_save.detrend(type='linear')
        data_array.append(tr1_save)
    if cor_events[r][9]==1:
        tr2a=stream2[int(cor_events[r][1])]
        dstart=UTCDateTime(cor_events[r][0])
        tr2_save = tr2a.slice(starttime=dstart, endtime=dstart+80)
        tr2_save.detrend(type='demean')
        tr2_save.detrend(type='linear')
        data_array.append(tr2_save)    
    if cor_events[r][10]==1:
        tr3a=stream3[int(cor_events[r][1])]
        dstart=UTCDateTime(cor_events[r][0])
        tr3_save = tr3a.slice(starttime=dstart, endtime=dstart+80)
        tr3_save.detrend(type='demean')
        tr3_save.detrend(type='linear')
        data_array.append(tr3_save)         
    if cor_events[r][11]==1:
        tr4a=stream4[int(cor_events[r][1])]
        dstart=UTCDateTime(cor_events[r][0])
        tr4_save = tr4a.slice(starttime=dstart, endtime=dstart+80)
        tr4_save.detrend(type='demean')
        tr4_save.detrend(type='linear')
        data_array.append(tr4_save)    
    if cor_events[r][12]==1:
        tr5a=stream5[int(cor_events[r][1])]
        dstart=UTCDateTime(cor_events[r][0])
        tr5_save = tr5a.slice(starttime=dstart, endtime=dstart+80)
        tr5_save.detrend(type='demean')
        tr5_save.detrend(type='linear')
        data_array.append(tr5_save)
    if cor_events[r][13]==1:
        tr6a=stream6[int(cor_events[r][1])] 
        dstart=UTCDateTime(cor_events[r][0])
        tr6_save = tr6a.slice(starttime=dstart, endtime=dstart+80)
        tr6_save.detrend(type='demean')
        tr6_save.detrend(type='linear')
        data_array.append(tr6_save)
        
        
#%%
print("End of Scan")
T2=time.clock()
elapsed= T2-T1
print('Time taken:', elapsed)
file_exp.close() 
file_cor2.close()
file_lb01.close()  
file_lb02.close()  
file_lb03.close()  
file_lb04.close()  
file_lb05.close()  
file_lb06.close()           
    
   
    
    
    
    
    
    
    
    