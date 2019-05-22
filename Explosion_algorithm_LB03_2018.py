#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 09:45:44 2019

@author: william
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 15:19:14 2018

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
from obspy import Stream

#T1= time.clock()

#%%
file_exp = open('/Users/william/Documents/scanner/all_stations/All_exp_scan_LB03_2018.txt','w')
file_cor2 = open('/Users/william/Documents/scanner/all_stations/All_correlated2_LB03_2018.txt','w')
#file_lb01= open('/Users/william/Documents/scanner/LB01/LB01_all_exp_scan_2_month_31_32.txt','w')
#file_lb02= open('/Users/william/Documents/scanner/LB02/LB02_all_exp_scan_2_month_31_32.txt','w')
file_lb03= open('/Users/william/Documents/scanner/LB03_2018.txt','w')
#file_lb04= open('/Users/william/Documents/scanner/LB04/LB04_all_exp_scan_2_month_31_32.txt','w')
#file_lb05= open('/Users/william/Documents/scanner/LB05/LB05_all_exp_scan_2_month_31_32.txt','w')
#file_lb06= open('/Users/william/Documents/scanner/LB06/LB06_all_exp_scan_2_month_31_32.txt','w')
#file_ls01= open('/Users/william/Documents/scanner/LS01/LS01_all_exp_scan_2_month_31_32.txt','w')
#file_ls02= open('/Users/william/Documents/scanner/LS02/LS02_all_exp_scan_2_month_31_32.txt','w')
#file_ls03= open('/Users/william/Documents/scanner/LS03/LS03_all_exp_scan_2_month_31_32.txt','w')
#file_ls04= open('/Users/william/Documents/scanner/LS04/LS04_all_exp_scan_2_month_31_32.txt','w')
#file_ls05= open('/Users/william/Documents/scanner/LS05/LS05_all_exp_scan_2_month_31_32.txt','w')
#file_ls06= open('/Users/william/Documents/scanner/LS06/LS06_all_exp_scan_2_month_31_32.txt','w')


   #%% constants  
step=4
shift=50*step
crite=0.6           # 0.7 in version 1
fmin=0.1
fmax=10           
#%% Reference stacked waveforms     

trs_e1,trs_e2,trs_e3,trs_e4,trs_e5,trs_e6=stacked_LBz_exp(fmin,fmax) 
#trs_s1,trs_s2,trs_s3,trs_s4,trs_s5,trs_s6=stacked_LS_exp(fmin,fmax) 
print('stacked waveform input')
#%% open empty streams for each station

#stream1 = Stream()
#stream2 = Stream()
stream3 = Stream()
#stream4 = Stream()
#stream5 = Stream()
#stream6 = Stream()
#stream1s = Stream()
#stream2s = Stream()
#stream3s = Stream()
#stream4s = Stream()
#stream5s = Stream()
#stream6s = Stream()


#%% Read in Each station
num_active=1
lb_num_active=1

stream3 = read("/Volumes/Wills_Disc/2018_seismic_data/STG8/HHZ.D/GI.STG8..HHZ.D.2018.*")

    
print('files read in') 

#%% Open matracies

Events_per_day=np.zeros(shape=(1,2))
Event_store=[]
All_events=[]


#%% Loop over whole time one day at a time
cor_events=np.zeros(shape=(1,20))
cor_count=0
for p in range(0,len(stream3)):
    file_exp.write("\n")
    file_cor2.write("\n")
    print('Day',p+1,'of',len(stream3))
    Events_today=np.zeros(shape=(1,14))
    todays_event=[]  
    correlated_today=[]                                    
    
    Number_of_events_today=0
    find_count=0    
    

    trace3=stream3[p]

        
# 
 #%% Station 3 scan  
    if trace3.stats.npts > 201:
        trace3.detrend(type='linear')
        trace3.detrend(type='demean')
        trace3.filter("bandpass", freqmin=fmin,freqmax=fmax)
        
        av_amp3=200
        
        begin = trace3.stats.starttime
        window=int(len(trs_e3)/100)
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
            if trc3.stats.npts < 100:
                break
            trc3.detrend(type='linear')
            trc3.detrend(type='demean')
            # only look for things above noise level 
            amp=abs(trc3.max())
            
            if amp > av_amp3:       
                peak3,cf3,bwid3,bwid253 = freq_info(trc3.data,trc3.stats.starttime,trc3.stats.endtime)      
                if (0.9 < cf3 < 4.25) and (0.2 < peak3 < 2.5) and (0.2 < bwid3 < 2.75):
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
  
                            row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,p])
                            
                            
                            todays_event.append(expt)
                            Event_store.append(row)
                            Events_today[find_count][0]=rt
                            Events_today[find_count][1]=p
                            Events_today[find_count][4]=1
                            Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
                            All_events.append(expt)
                            find_count += 1
                        
                        else:
                            Events_today[ind][4]=1

                        
                                   
                        
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
                        cor_count += 1
                        Number_of_events_today += 1
                        lastcorr=0    
                        add += window-step  +30             #skip 60s to avoid double catch

            else:
                add += window-step
                lastcorr=0 
    else:
        switch3=0
                
                 #%% Station 4 scan  
#    if trace4.stats.npts > 201:
#        trace4.detrend(type='linear')
#        trace4.detrend(type='demean')
#        trace4.filter("bandpass", freqmin=fmin,freqmax=fmax)
#        
#        av_amp4=450000
#
#    
#        begin = trace4.stats.starttime
#        window=int(len(trs_e4)/100)
#        add=0
#        lastcorr=0  
#        expt=0
#        dl=trace4.stats.endtime - trace4.stats.starttime     
#        for t in range(0,int(trace4.stats.npts/100), step):
#            view_start = begin + t + add
#            view_end = view_start + window
#            if view_start > trace4.stats.endtime - window :
#                break
#    
#            trc4 = trace4.slice(starttime = view_start  , endtime= view_end)
#            if trc4.stats.npts < 100:
#                break
#            trc4.detrend(type='linear')
#            trc4.detrend(type='demean')
#            # only look for things above noise level 
#            amp=abs(trc4.max())
#            
#            if amp > av_amp4:       
#                peak4,cf4,bwid4,bwid254 = freq_info(trc4.data,trc4.stats.starttime,trc4.stats.endtime)      
#                if (0.9 < cf4 < 5) and (0.4 < peak4 < 3.5) and (0.3 < bwid4 < 4): 
#                         #correlate between data and earthquake envelopes 
#                    trc4_e = obspy.signal.filter.envelope(trc4.data) 
#                    top_v,top,corell = corel(trs_e4,trc4_e,shift)
#                   
#                    if 0.0 in corell : # in the case of missing data - which will cause function to crash
#                        break            
# 
#                    expt_new = view_start - (top/100) +10
#
#                    if abs(top_v) > crite and top_v > lastcorr:
#                            lastcorr = top_v
#                            expt = expt_new                                    
#                        
#                    if lastcorr > top_v:
#
#                        file_lb04.write(str(expt))
#                        file_lb04.write("\n")
#                        
#                        rt=expt.timestamp	
#                        jd=expt.julday
#                        yr=expt.year
#                        mo=expt.month
#                        da=expt.day
#                        hr=expt.hour
#                        mi=expt.minute
#                        se=expt.second
#                        ms=expt.microsecond/10000
#                        
#                        near,ind=find_nearest(Events_today[:,0], rt)
#                        if abs(rt-near) > 60: #60s
#  
#                            row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,x])
#                            
#                            todays_event.append(expt)
#                            Event_store.append(row)
#                            Events_today[find_count][0]=rt
#                            Events_today[find_count][1]=p
#                            Events_today[find_count][5]=1
#                            Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            All_events.append(expt)
#                            find_count += 1
#                        else:
#                            Events_today[ind][5]=1
#                        
#                        if lb_num_active[p] > 1:
#                            if abs(rt-near) <= 20:
#                                near_cor,inde=find_nearest(cor_events[:,0], rt)
#                                if abs(rt-near_cor) >  60:
#                                    cor_events[cor_count][0]=rt
#                                    cor_events[cor_count][1]="%1d" % p
#                                    cor_events[cor_count][2]="%1d" % yr
#                                    cor_events[cor_count][3]="%1d" % jd
#                                    cor_events[cor_count][4]="%1d" % hr
#                                    cor_events[cor_count][5]="%1d" % mi
#                                    cor_events[cor_count][6]="%1d" % se
#                                    cor_events[cor_count][7]="%1d" % ms
#                                    cor_events[cor_count][8]="%1d" % Events_today[ind][2]
#                                    cor_events[cor_count][9]="%1d" % Events_today[ind][3]
#                                    cor_events[cor_count][10]="%1d" % Events_today[ind][4]
#                                    cor_events[cor_count][11]="%1d" % Events_today[ind][5]
#                                    cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                    correlated_today.append(expt)
#                                    cor_count += 1
#                                    Number_of_events_today += 1
#                                else:
#                                    cor_events[inde][11]="%1d" % 1
#                        else:
#                            cor_events[cor_count][0]=rt
#                            cor_events[cor_count][1]="%1d" % p
#                            cor_events[cor_count][2]="%1d" % yr
#                            cor_events[cor_count][3]="%1d" % jd
#                            cor_events[cor_count][4]="%1d" % hr
#                            cor_events[cor_count][5]="%1d" % mi
#                            cor_events[cor_count][6]="%1d" % se
#                            cor_events[cor_count][7]="%1d" % ms
#                            cor_events[cor_count][11]="%1d" % 1
#                            cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            correlated_today.append(expt)
#                            cor_count += 1
#                            Number_of_events_today += 1
#                            
#                        lastcorr=0    
#                        add += window-step  +30           #skip 60s to avoid double catch
#
#            else:
#                add += window-step
#                lastcorr=0 
#    else:
#        switch4=0
#                
#    #%% Station 5 scan  
#    if trace5.stats.npts > 201:
#        trace5.detrend(type='linear')
#        trace5.detrend(type='demean')
#        trace5.filter("bandpass", freqmin=fmin,freqmax=fmax)
#        
#        av_amp5=300000
#        
#        begin = trace5.stats.starttime
#        window=int(len(trs_e5)/100)
#        add=0
#        lastcorr=0 
#        expt=0
#        dl=trace5.stats.endtime - trace5.stats.starttime    
#        for t in range(0,int(trace5.stats.npts/100), step):
#            view_start = begin + t + add
#            view_end = view_start + window
#            if view_start > trace5.stats.endtime - window :
#                break
#    
#            trc5 = trace5.slice(starttime = view_start  , endtime= view_end)
#            if trc5.stats.npts < 100:
#                break
#            trc5.detrend(type='linear')
#            trc5.detrend(type='demean')
#            # only look for things above noise level 
#            amp=abs(trc5.max())
#            
#            if amp > av_amp5:       
#                peak5,cf5,bwid5,bwid255 = freq_info(trc5.data,trc5.stats.starttime,trc5.stats.endtime)      
#                if (0.9 < cf5 < 5) and (0.4 < peak5 < 3.5) and (0.3 < bwid5 < 3.5): 
#                         #correlate between data and earthquake envelopes 
#                    trc5_e = obspy.signal.filter.envelope(trc5.data) 
#                    top_v,top,corell = corel(trs_e5,trc5_e,shift)
#                    
#                    if 0.0 in corell : # in the case of missing data - which will cause function to crash
#                        break            
# 
#                    expt_new = view_start - (top/100) +10
#
#                    if abs(top_v) > crite and top_v > lastcorr:
#                            lastcorr = top_v
#                            expt = expt_new                                    
#                        
#                    if lastcorr > top_v:
#
#                        file_lb05.write(str(expt))
#                        file_lb05.write("\n")
#                        
#                        rt=expt.timestamp	
#                        jd=expt.julday
#                        yr=expt.year
#                        mo=expt.month
#                        da=expt.day
#                        hr=expt.hour
#                        mi=expt.minute
#                        se=expt.second
#                        ms=expt.microsecond/10000
#                        near,ind=find_nearest(Events_today[:,0], rt)
#                        
#                        if abs(rt-near) > 60: #60s
#  
#                            row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,x])
#                            
#                            todays_event.append(expt)
#                            Event_store.append(row)
#                            Events_today[find_count][0]=rt
#                            Events_today[find_count][1]=p
#                            Events_today[find_count][6]=1
#                            Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            All_events.append(expt)
#                            find_count += 1
#                        else:
#                            Events_today[ind][6]=1
#                        
#                        if lb_num_active[p] > 1:
#                            if abs(rt-near) <= 20:
#                                near_cor,inde=find_nearest(cor_events[:,0], rt)
#                                if abs(rt-near_cor) >  60:
#                                    cor_events[cor_count][0]=rt
#                                    cor_events[cor_count][1]="%1d" % p
#                                    cor_events[cor_count][2]="%1d" % yr
#                                    cor_events[cor_count][3]="%1d" % jd
#                                    cor_events[cor_count][4]="%1d" % hr
#                                    cor_events[cor_count][5]="%1d" % mi
#                                    cor_events[cor_count][6]="%1d" % se
#                                    cor_events[cor_count][7]="%1d" % ms
#                                    cor_events[cor_count][8]="%1d" % Events_today[ind][2]
#                                    cor_events[cor_count][9]="%1d" % Events_today[ind][3]
#                                    cor_events[cor_count][10]="%1d" % Events_today[ind][4]
#                                    cor_events[cor_count][11]="%1d" % Events_today[ind][5]
#                                    cor_events[cor_count][12]="%1d" % Events_today[ind][6]
#                                    cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                    correlated_today.append(expt)
#                                    cor_count += 1
#                                    Number_of_events_today += 1
#                                else:
#                                    cor_events[inde][12]="%1d" % 1
#
#                        else:
#                            cor_events[cor_count][0]= rt
#                            cor_events[cor_count][1]="%1d" % p
#                            cor_events[cor_count][2]="%1d" % yr
#                            cor_events[cor_count][3]="%1d" % jd
#                            cor_events[cor_count][4]="%1d" % hr
#                            cor_events[cor_count][5]="%1d" % mi
#                            cor_events[cor_count][6]="%1d" % se
#                            cor_events[cor_count][7]="%1d" % ms
#                            cor_events[cor_count][12]="%1d" % 1
#                            cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            correlated_today.append(expt)
#                            cor_count += 1
#                            Number_of_events_today += 1
#                        
#                        lastcorr=0
#                        add += window-step  +30           #skip 60s to avoid double catch
#
#            else:
#                add += window-step
#                lastcorr=0
#    else:
#        switch5=0
#                
#    #%% Station 6 scan  
#    if trace6.stats.npts > 201:
#        trace6.detrend(type='linear')
#        trace6.detrend(type='demean')
#        trace6.filter("bandpass", freqmin=fmin,freqmax=fmax)
#        
#        av_amp6=700
#    
#        begin = trace6.stats.starttime
#        window=int(len(trs_e6)/100)
#        add=0
#        lastcorr=0 
#        expt=0
#        dl=trace6.stats.endtime - trace6.stats.starttime    
#        for t in range(0,int(trace6.stats.npts/100), step):
#            view_start = begin + t + add
#            view_end = view_start + window
#            if view_start > trace6.stats.endtime - window :
#                break
#    
#            trc6 = trace6.slice(starttime = view_start  , endtime= view_end)
#            if trc6.stats.npts < 100:
#                break
#            trc6.detrend(type='linear')
#            trc6.detrend(type='demean')
#            # only look for things above noise level 
#            amp=abs(trc6.max())
#            
#            if amp > av_amp6:       
#                peak6,cf6,bwid6,bwid256 = freq_info(trc6.data,trc6.stats.starttime,trc6.stats.endtime)      
#                if (0.65 < cf6 < 2.75) and (0.2 < peak6 < 2.5) and (0.3 < bwid6 < 3): 
#                         #correlate between data and earthquake envelopes 
#                    trc6_e = obspy.signal.filter.envelope(trc6.data) 
#                    top_v,top,corell = corel(trs_e6,trc6_e,shift)
#                    
#                    if 0.0 in corell : # in the case of missing data - which will cause function to crash
#                        break            
# 
#                    expt_new = view_start - (top/100) +10
#
#                    if abs(top_v) > crite and top_v > lastcorr:
#                            lastcorr = top_v
#                            expt = expt_new                                    
#                        
#                    if lastcorr > top_v:
#
#                        file_lb06.write(str(expt))
#                        file_lb06.write("\n")
#                        
#                        rt=expt.timestamp	
#                        jd=expt.julday
#                        yr=expt.year
#                        mo=expt.month
#                        da=expt.day
#                        hr=expt.hour
#                        mi=expt.minute
#                        se=expt.second
#                        ms=expt.microsecond/10000
#                        
#                        near,ind=find_nearest(Events_today[:,0], rt)
#                        if abs(rt-near) > 60: #60s
# 
#                            row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl,x])
#                            
#                            todays_event.append(expt)
#                            Event_store.append(row)
#                            Events_today[find_count][0]=rt
#                            Events_today[find_count][1]=p
#                            Events_today[find_count][7]=1
#                            Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            All_events.append(expt)
#                            find_count += 1
#                        else:
#                            Events_today[ind][7]=1
#                        
#                        if lb_num_active[p] > 1:
#                            if abs(rt-near) <= 20:
#                                near_cor,inde=find_nearest(cor_events[:,0], rt)
#                                if abs(rt-near_cor) >  60:
#                                    cor_events[cor_count][0]=rt
#                                    cor_events[cor_count][1]="%1d" % p
#                                    cor_events[cor_count][2]="%1d" % yr
#                                    cor_events[cor_count][3]="%1d" % jd
#                                    cor_events[cor_count][4]="%1d" % hr
#                                    cor_events[cor_count][5]="%1d" % mi
#                                    cor_events[cor_count][6]="%1d" % se
#                                    cor_events[cor_count][7]="%1d" % ms
#                                    cor_events[cor_count][8]="%1d" % Events_today[ind][2]
#                                    cor_events[cor_count][9]="%1d" % Events_today[ind][3]
#                                    cor_events[cor_count][10]="%1d" % Events_today[ind][4]
#                                    cor_events[cor_count][11]="%1d" % Events_today[ind][5]
#                                    cor_events[cor_count][12]="%1d" % Events_today[ind][6]
#                                    cor_events[cor_count][13]="%1d" % Events_today[ind][7]
#                                    cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                    correlated_today.append(expt)
#                                    cor_count += 1
#                                    Number_of_events_today += 1
#                                else:
#                                    cor_events[inde][13]="%1d" % 1
#                        else:
#                            cor_events[cor_count][0]=rt
#                            cor_events[cor_count][1]="%1d" % p
#                            cor_events[cor_count][2]="%1d" % yr
#                            cor_events[cor_count][3]="%1d" % jd
#                            cor_events[cor_count][4]="%1d" % hr
#                            cor_events[cor_count][5]="%1d" % mi
#                            cor_events[cor_count][6]="%1d" % se
#                            cor_events[cor_count][7]="%1d" % ms
#                            cor_events[cor_count][13]="%1d" % 1
#                            cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            correlated_today.append(expt)                           
#                            cor_count += 1
#                            Number_of_events_today += 1
#                        lastcorr=0                                        
#                        add += window-step  +30           #skip 60s to avoid double catch
#
#            else:
#                add += window-step
#                lastcorr=0 
#    else:
#        switch6=0
#
##%% Station LS01 scan    
#    if trace1s.stats.npts > 201: 
#        trace1s.detrend(type='linear')
#        trace1s.detrend(type='demean')
#        trace1s.filter("bandpass", freqmin=fmin,freqmax=fmax)
#        
#        av_amp1s=5000
#    
#        begin = trace1s.stats.starttime
#        window=int(len(trs_s1)/100)
#        add=0
#        lastcorr=0
#        expt=0
#        dl=trace1s.stats.endtime - trace1s.stats.starttime 
#        # Loop over single day, scanning for waveforms
#        for t in range(0,int(trace1s.stats.npts/100), step):
#            view_start = begin + t + add
#            view_end = view_start + window
#            if view_start > trace1s.stats.endtime - window :
#                break
#    
#            trcs = trace1s.slice(starttime = view_start  , endtime= view_end)
#            if trcs.stats.npts < 100:
#                break
#            trcs.detrend(type='linear')
#            trcs.detrend(type='demean')
#            # only look for things above noise level 
#            amp=abs(trcs.max())
#            
#            if amp > av_amp1s:                
#                peak,cf,bwid,bwid25 = freq_info(trcs.data,trcs.stats.starttime,trcs.stats.endtime)         
#                if (1.5 < cf < 3.75) and (1.1 < peak < 3.25) and (0.1 < bwid < 2.75): 
#                         #correlate between data and earthquake envelopes 
#                    trcs_e = obspy.signal.filter.envelope(trcs.data) 
#                    top_v,top,corell = corel(trs_s1,trcs_e,shift)
#                    
#                    if 0.0 in corell : # in the case of missing data - which will cause function to crash
#                        break   
#                    
#                    expt_new = view_start - (top/100) +10
#
#                    if abs(top_v) > crite and top_v > lastcorr:
#                            lastcorr = top_v
#                            expt = expt_new                                    
#                        
#                    if lastcorr > top_v:
#
#                        file_ls01.write(str(expt))
#                        file_ls01.write("\n")
#                        
#                        rt=expt.timestamp	
#                        jd=expt.julday
#                        yr=expt.year
#                        mo=expt.month
#                        da=expt.day
#                        hr=expt.hour
#                        mi=expt.minute
#                        se=expt.second
#                        ms=expt.microsecond/10000
#                        near,ind=find_nearest(Events_today[:,0], rt)
#
#                        if num_active[p] > 1:
#                            if (lb_num_active[p] == 0) and (abs(rt-near) > 60):
#                                Events_today[find_count][0]=rt
#                                Events_today[find_count][1]=p
#                                Events_today[find_count][8]=1
#                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                All_events.append(expt)
#                                find_count += 1
#                            if abs(rt-near) <= 20:
#                                Events_today[ind][8]=1
#                                if (lb_num_active[p] == 0) or ((lb_num_active[p] > 0) and ((Events_today[ind][2] == 1) or (Events_today[ind][4] == 1))):
#                                    near_cor,inde=find_nearest(cor_events[:,0], rt)
#                                    if abs(rt-near_cor) >  60:
#                                        cor_events[cor_count][0]=rt
#                                        cor_events[cor_count][1]="%1d" % p
#                                        cor_events[cor_count][2]="%1d" % yr
#                                        cor_events[cor_count][3]="%1d" % jd
#                                        cor_events[cor_count][4]="%1d" % hr
#                                        cor_events[cor_count][5]="%1d" % mi
#                                        cor_events[cor_count][6]="%1d" % se
#                                        cor_events[cor_count][7]="%1d" % ms
#                                        cor_events[cor_count][8]="%1d" % Events_today[ind][2]
#                                        cor_events[cor_count][9]="%1d" % Events_today[ind][3]
#                                        cor_events[cor_count][10]="%1d" % Events_today[ind][4]
#                                        cor_events[cor_count][11]="%1d" % Events_today[ind][5]
#                                        cor_events[cor_count][12]="%1d" % Events_today[ind][6]
#                                        cor_events[cor_count][13]="%1d" % Events_today[ind][7]
#                                        cor_events[cor_count][14]="%1d" % Events_today[ind][8]
#                                        cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        correlated_today.append(expt)
#                                        cor_count += 1
#                                        Number_of_events_today += 1
#                                    else:
#                                        cor_events[inde][14]="%1d" % 1
#                                        
#                                   
#                        else:
#                            cor_events[cor_count][0]=rt
#                            cor_events[cor_count][1]="%1d" % p
#                            cor_events[cor_count][2]="%1d" % yr
#                            cor_events[cor_count][3]="%1d" % jd
#                            cor_events[cor_count][4]="%1d" % hr
#                            cor_events[cor_count][5]="%1d" % mi
#                            cor_events[cor_count][6]="%1d" % se
#                            cor_events[cor_count][7]="%1d" % ms
#                            cor_events[cor_count][14]="%1d" % 1
#                            cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            correlated_today.append(expt)
#                            cor_count += 1
#                            Number_of_events_today += 1
#                        lastcorr=0    
#                        add += window-step  +30             #skip 60s to avoid double catch
#
#            else:
#                add += window-step
#                lastcorr=0 
#    else:
#        switch1s=0
# #%% Station 2 scan  
#    if trace2s.stats.npts > 201:
#        trace2s.detrend(type='linear')
#        trace2s.detrend(type='demean')
#        trace2s.filter("bandpass", freqmin=fmin,freqmax=fmax)
#        
#        av_amp2s=1200
#        
#        begin = trace2s.stats.starttime
#        window=int(len(trs_s2)/100)
#        add=0
#        lastcorr=0     
#        expt=0
#        dl=trace2s.stats.endtime - trace2s.stats.starttime     
#        for t in range(0,int(trace2s.stats.npts/100), step):
#            view_start = begin + t + add
#            view_end = view_start + window
#            if view_start > trace2s.stats.endtime - window :
#                break
#    
#            trc2s = trace2s.slice(starttime = view_start  , endtime= view_end)
#            if trc2s.stats.npts < 100:
#                break
#            trc2s.detrend(type='linear')
#            trc2s.detrend(type='demean')
#            # only look for things above noise level 
#            amp=abs(trc2s.max())
#            
#            if amp > av_amp2s:  
#                peak2,cf2,bwid2,bwid252 = freq_info(trc2s.data,trc2s.stats.starttime,trc2s.stats.endtime)      
#                if (1.5 < cf2 < 4) and (1.1 < peak2 < 2.75) and (0.1 < bwid2 < 2.75):
#                         #correlate between data and earthquake envelopes 
#                    trc2s_e = obspy.signal.filter.envelope(trc2s.data) 
#                    top_v,top,corell = corel(trs_s2,trc2s_e,shift)
#                    
#                    if 0.0 in corell : # in the case of missing data - which will cause function to crash
#                        break            
#                    
#                    expt_new = view_start - (top/100) +10
#        #            if the correlation is positive, exp found
#                    if abs(top_v) > crite and top_v > lastcorr:
#                        lastcorr = top_v
#                        expt = expt_new  
#
#                    if lastcorr > top_v:
#                    
#                        file_ls02.write(str(expt))
#                        file_ls02.write("\n")
##                                
#                        rt=expt.timestamp	
#                        jd=expt.julday
#                        yr=expt.year
#                        mo=expt.month
#                        da=expt.day
#                        hr=expt.hour
#                        mi=expt.minute
#                        se=expt.second
#                        ms=expt.microsecond/10000
#                        near,ind=find_nearest(Events_today[:,0], rt)
#
#                        if num_active[p] > 1:
#                            if (lb_num_active[p] == 0) and (abs(rt-near) > 60):
#                                Events_today[find_count][0]=rt
#                                Events_today[find_count][1]=p
#                                Events_today[find_count][9]=1
#                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                All_events.append(expt)
#                                find_count += 1
#                            if abs(rt-near) <= 20:
#                                Events_today[ind][9]=1
#                                if (lb_num_active[p] == 0) or ((lb_num_active[p] > 0) and ((Events_today[ind][2] == 1) or (Events_today[ind][4] == 1))):
#                                    near_cor,inde=find_nearest(cor_events[:,0], rt)
#                                    if abs(rt-near_cor) >  60:
#                                        cor_events[cor_count][0]=rt
#                                        cor_events[cor_count][1]="%1d" % p
#                                        cor_events[cor_count][2]="%1d" % yr
#                                        cor_events[cor_count][3]="%1d" % jd
#                                        cor_events[cor_count][4]="%1d" % hr
#                                        cor_events[cor_count][5]="%1d" % mi
#                                        cor_events[cor_count][6]="%1d" % se
#                                        cor_events[cor_count][7]="%1d" % ms
#                                        cor_events[cor_count][8]="%1d" % Events_today[ind][2]
#                                        cor_events[cor_count][9]="%1d" % Events_today[ind][3]
#                                        cor_events[cor_count][10]="%1d" % Events_today[ind][4]
#                                        cor_events[cor_count][11]="%1d" % Events_today[ind][5]
#                                        cor_events[cor_count][12]="%1d" % Events_today[ind][6]
#                                        cor_events[cor_count][13]="%1d" % Events_today[ind][7]
#                                        cor_events[cor_count][14]="%1d" % Events_today[ind][8]
#                                        cor_events[cor_count][15]="%1d" % Events_today[ind][9]
#                                        cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        correlated_today.append(expt)
#                                        cor_count += 1
#                                        Number_of_events_today += 1
#                                    else:
#                                        cor_events[inde][15]="%1d" % 1
#                                   
#                        else:
#                            cor_events[cor_count][0]=rt
#                            cor_events[cor_count][1]="%1d" % p
#                            cor_events[cor_count][2]="%1d" % yr
#                            cor_events[cor_count][3]="%1d" % jd
#                            cor_events[cor_count][4]="%1d" % hr
#                            cor_events[cor_count][5]="%1d" % mi
#                            cor_events[cor_count][6]="%1d" % se
#                            cor_events[cor_count][7]="%1d" % ms
#                            cor_events[cor_count][15]="%1d" % 1
#                            cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            correlated_today.append(expt)
#                            cor_count += 1
#                            Number_of_events_today += 1
#                        lastcorr=0    
#                        add += window-step +30              #skip 60s to avoid double catch
#
#            else:
#                add += window-step
#                lastcorr=0 
#    else:
#        switch2s=0 
# #%% Station 3 scan  
#    if trace3s.stats.npts > 201:
#        trace3s.detrend(type='linear')
#        trace3s.detrend(type='demean')
#        trace3s.filter("bandpass", freqmin=fmin,freqmax=fmax)
#        
#        av_amp3s=2000
#    
#        begin = trace3s.stats.starttime
#        window=int(len(trs_s3)/100)
#        add=0
#        lastcorr=0  
#        expt=0
#        dl=trace3s.stats.endtime - trace3s.stats.starttime     
#        for t in range(0,int(trace3s.stats.npts/100), step):
#            view_start = begin + t + add
#            view_end = view_start + window
#            if view_start > trace3s.stats.endtime - window :
#                break
#    
#            trc3s = trace3s.slice(starttime = view_start  , endtime= view_end)
#            if trc3s.stats.npts < 100:
#                break
#            trc3s.detrend(type='linear')
#            trc3s.detrend(type='demean')
#            # only look for things above noise level 
#            amp=abs(trc3s.max())
#            
#            if amp > av_amp3s:       
#                peak3,cf3,bwid3,bwid253 = freq_info(trc3s.data,trc3s.stats.starttime,trc3s.stats.endtime)      
#                if (1.75 < cf3 < 5.5) and (1.1 < peak3 < 3.5) and (0.1 < bwid3 < 3.75):
#                         #correlate between data and earthquake envelopes 
#                    trc3s_e = obspy.signal.filter.envelope(trc3s.data) 
#                    top_v,top,corell = corel(trs_s3,trc3s_e,shift)
#                    
#                    if 0.0 in corell : # in the case of missing data - which will cause function to crash
#                        break            
# 
#                    expt_new = view_start - (top/100) +10
#
#                    if abs(top_v) > crite and top_v > lastcorr:
#                            lastcorr = top_v
#                            expt = expt_new                                    
#                        
#                    if lastcorr > top_v:
#                        
#                        file_ls03.write(str(expt))
#                        file_ls03.write("\n")
#                        
#                        rt=expt.timestamp	
#                        jd=expt.julday
#                        yr=expt.year
#                        mo=expt.month
#                        da=expt.day
#                        hr=expt.hour
#                        mi=expt.minute
#                        se=expt.second
#                        ms=expt.microsecond/10000
#                        near,ind=find_nearest(Events_today[:,0], rt)
#
#                        if num_active[p] > 1:
#                            if (lb_num_active[p] == 0) and (abs(rt-near) > 60):
#                                Events_today[find_count][0]=rt
#                                Events_today[find_count][1]=p
#                                Events_today[find_count][10]=1
#                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                All_events.append(expt)
#                                find_count += 1
#                            if abs(rt-near) <= 20:
#                                Events_today[ind][10]=1
#                                if (lb_num_active[p] == 0) or ((lb_num_active[p] > 0) and ((Events_today[ind][2] == 1) or (Events_today[ind][4] == 1))):
#                                    near_cor,inde=find_nearest(cor_events[:,0], rt)
#                                    if abs(rt-near_cor) >  60:
#                                        cor_events[cor_count][0]=rt
#                                        cor_events[cor_count][1]="%1d" % p
#                                        cor_events[cor_count][2]="%1d" % yr
#                                        cor_events[cor_count][3]="%1d" % jd
#                                        cor_events[cor_count][4]="%1d" % hr
#                                        cor_events[cor_count][5]="%1d" % mi
#                                        cor_events[cor_count][6]="%1d" % se
#                                        cor_events[cor_count][7]="%1d" % ms
#                                        cor_events[cor_count][8]="%1d" % Events_today[ind][2]
#                                        cor_events[cor_count][9]="%1d" % Events_today[ind][3]
#                                        cor_events[cor_count][10]="%1d" % Events_today[ind][4]
#                                        cor_events[cor_count][11]="%1d" % Events_today[ind][5]
#                                        cor_events[cor_count][12]="%1d" % Events_today[ind][6]
#                                        cor_events[cor_count][13]="%1d" % Events_today[ind][7]
#                                        cor_events[cor_count][14]="%1d" % Events_today[ind][8]
#                                        cor_events[cor_count][15]="%1d" % Events_today[ind][9]
#                                        cor_events[cor_count][16]="%1d" % Events_today[ind][10]
#                                        cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        correlated_today.append(expt)
#                                        cor_count += 1
#                                        Number_of_events_today += 1
#                                    else:
#                                        cor_events[inde][16]="%1d" % 1
#                                   
#                        else:
#                            cor_events[cor_count][0]=rt
#                            cor_events[cor_count][1]="%1d" % p
#                            cor_events[cor_count][2]="%1d" % yr
#                            cor_events[cor_count][3]="%1d" % jd
#                            cor_events[cor_count][4]="%1d" % hr
#                            cor_events[cor_count][5]="%1d" % mi
#                            cor_events[cor_count][6]="%1d" % se
#                            cor_events[cor_count][7]="%1d" % ms
#                            cor_events[cor_count][16]="%1d" % 1
#                            cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            correlated_today.append(expt)
#                            cor_count += 1
#                            Number_of_events_today += 1
#                        lastcorr=0    
#                        add += window-step +30           #skip 60s to avoid double catch
#
#            else:
#                add += window-step
#                lastcorr=0 
#    else:
#        switch3s=0
#                
#                 #%% Station 4 scan  
#    if trace4s.stats.npts > 201:
#        trace4s.detrend(type='linear')
#        trace4s.detrend(type='demean')
#        trace4s.filter("bandpass", freqmin=fmin,freqmax=fmax)
#        if trace4s.stats.starttime.timestamp < 1448841601.0:
#            av_amp4s=12000
#        else:
#            av_amp4s=600
#        
#        begin = trace4s.stats.starttime
#        window=int(len(trs_s4)/100)
#        add=0
#        lastcorr=0  
#        expt=0
#        dl=trace4s.stats.endtime - trace4s.stats.starttime     
#        for t in range(0,int(trace4s.stats.npts/100), step):
#            view_start = begin + t + add
#            view_end = view_start + window
#            if view_start > trace4s.stats.endtime - window :
#                break
#    
#            trc4s = trace4s.slice(starttime = view_start  , endtime= view_end)
#            if trc4s.stats.npts < 100:
#                break
#            trc4s.detrend(type='linear')
#            trc4s.detrend(type='demean')
#            # only look for things above noise level 
#            amp=abs(trc4s.max())
#            
#            if amp > av_amp4s:       
#                peak4,cf4,bwid4,bwid254 = freq_info(trc4s.data,trc4s.stats.starttime,trc4s.stats.endtime)      
#                if (2.25 < cf4 < 5) and (1.1 < peak4 < 4.75) and (0.1 < bwid4 < 4.75): 
#                         #correlate between data and earthquake envelopes 
#                    trc4s_e = obspy.signal.filter.envelope(trc4s.data) 
#                    top_v,top,corell = corel(trs_s4,trc4s_e,shift)
#                   
#                    if 0.0 in corell : # in the case of missing data - which will cause function to crash
#                        break            
# 
#                    expt_new = view_start - (top/100) +10
#
#                    if abs(top_v) > crite and top_v > lastcorr:
#                            lastcorr = top_v
#                            expt = expt_new                                    
#                        
#                    if lastcorr > top_v:
#
#                        file_ls04.write(str(expt))
#                        file_ls04.write("\n")
#                        
#                        rt=expt.timestamp	
#                        jd=expt.julday
#                        yr=expt.year
#                        mo=expt.month
#                        da=expt.day
#                        hr=expt.hour
#                        mi=expt.minute
#                        se=expt.second
#                        ms=expt.microsecond/10000
#                        near,ind=find_nearest(Events_today[:,0], rt)
#
#                        if num_active[p] > 1:
#                            if (lb_num_active[p] == 0) and (abs(rt-near) > 60):
#                                Events_today[find_count][0]=rt
#                                Events_today[find_count][1]=p
#                                Events_today[find_count][11]=1
#                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                All_events.append(expt)
#                                find_count += 1
#                            if abs(rt-near) <= 20:
#                                Events_today[ind][11]=1
#                                if (lb_num_active[p] == 0) or ((lb_num_active[p] > 0) and ((Events_today[ind][2] == 1) or (Events_today[ind][4] == 1))):
#                                    near_cor,inde=find_nearest(cor_events[:,0], rt)
#                                    if abs(rt-near_cor) >  60:
#                                        cor_events[cor_count][0]=rt
#                                        cor_events[cor_count][1]="%1d" % p
#                                        cor_events[cor_count][2]="%1d" % yr
#                                        cor_events[cor_count][3]="%1d" % jd
#                                        cor_events[cor_count][4]="%1d" % hr
#                                        cor_events[cor_count][5]="%1d" % mi
#                                        cor_events[cor_count][6]="%1d" % se
#                                        cor_events[cor_count][7]="%1d" % ms
#                                        cor_events[cor_count][8]="%1d" % Events_today[ind][2]
#                                        cor_events[cor_count][9]="%1d" % Events_today[ind][3]
#                                        cor_events[cor_count][10]="%1d" % Events_today[ind][4]
#                                        cor_events[cor_count][11]="%1d" % Events_today[ind][5]
#                                        cor_events[cor_count][12]="%1d" % Events_today[ind][6]
#                                        cor_events[cor_count][13]="%1d" % Events_today[ind][7]
#                                        cor_events[cor_count][14]="%1d" % Events_today[ind][8]
#                                        cor_events[cor_count][15]="%1d" % Events_today[ind][9]
#                                        cor_events[cor_count][16]="%1d" % Events_today[ind][10]
#                                        cor_events[cor_count][17]="%1d" % Events_today[ind][11]
#                                        cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        correlated_today.append(expt)
#                                        cor_count += 1
#                                        Number_of_events_today += 1
#                                    else:
#                                        cor_events[inde][17]="%1d" % 1
#                        else:
#                            cor_events[cor_count][0]=rt
#                            cor_events[cor_count][1]="%1d" % p
#                            cor_events[cor_count][2]="%1d" % yr
#                            cor_events[cor_count][3]="%1d" % jd
#                            cor_events[cor_count][4]="%1d" % hr
#                            cor_events[cor_count][5]="%1d" % mi
#                            cor_events[cor_count][6]="%1d" % se
#                            cor_events[cor_count][7]="%1d" % ms
#                            cor_events[cor_count][17]="%1d" % 1
#                            cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            correlated_today.append(expt)
#                            cor_count += 1
#                            Number_of_events_today += 1
#                            
#                        lastcorr=0    
#                        add += window-step   +30           #skip 60s to avoid double catch
#
#            else:
#                add += window-step
#                lastcorr=0 
#    else:
#        switch4s=0
##                
#    #%% Station 5 scan  
#    if trace5s.stats.npts > 201:
#        trace5s.detrend(type='linear')
#        trace5s.detrend(type='demean')
#        trace5s.filter("bandpass", freqmin=fmin,freqmax=fmax)
#    
#        av_amp5s=2500
#        
#        begin = trace5s.stats.starttime
#        window=int(len(trs_s5)/100)
#        add=0
#        lastcorr=0 
#        expt=0
#        dl=trace5s.stats.endtime - trace5s.stats.starttime    
#        for t in range(0,int(trace5s.stats.npts/100), step):
#            view_start = begin + t + add
#            view_end = view_start + window
#            if view_start > trace5s.stats.endtime - window :
#                break
#    
#            trc5s = trace5s.slice(starttime = view_start  , endtime= view_end)
#            if trc5s.stats.npts < 100:
#                break
#            trc5s.detrend(type='linear')
#            trc5s.detrend(type='demean')
#            # only look for things above noise level 
#            amp=abs(trc5s.max())
#            
#            if amp > av_amp5s:       
#                peak5,cf5,bwid5,bwid255 = freq_info(trc5s.data,trc5s.stats.starttime,trc5s.stats.endtime)      
#                if (1.5 < cf5 < 5) and (1.1 < peak5 < 3) and (0.1 < bwid5 < 2.5): 
#                         #correlate between data and earthquake envelopes 
#                    trc5s_e = obspy.signal.filter.envelope(trc5s.data) 
#                    top_v,top,corell = corel(trs_s5,trc5s_e,shift)
#                    
#                    if 0.0 in corell : # in the case of missing data - which will cause function to crash
#                        break            
# 
#                    expt_new = view_start - (top/100) +10
#
#                    if abs(top_v) > crite and top_v > lastcorr:
#                            lastcorr = top_v
#                            expt = expt_new                                    
#                        
#                    if lastcorr > top_v:
#
#                        file_ls05.write(str(expt))
#                        file_ls05.write("\n")
#                        
#                        rt=expt.timestamp	
#                        jd=expt.julday
#                        yr=expt.year
#                        mo=expt.month
#                        da=expt.day
#                        hr=expt.hour
#                        mi=expt.minute
#                        se=expt.second
#                        ms=expt.microsecond/10000
#                        near,ind=find_nearest(Events_today[:,0], rt)
#
#                        if num_active[p] > 1:
#                            if (lb_num_active[p] == 0) and (abs(rt-near) > 60):
#                                Events_today[find_count][0]=rt
#                                Events_today[find_count][1]=p
#                                Events_today[find_count][12]=1
#                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                All_events.append(expt)
#                                find_count += 1
#                            if abs(rt-near) <= 20:
#                                Events_today[ind][12]=1
#                                if (lb_num_active[p] == 0) or ((lb_num_active[p] > 0) and ((Events_today[ind][2] == 1) or (Events_today[ind][4] == 1))):
#                                    near_cor,inde=find_nearest(cor_events[:,0], rt)
#                                    if abs(rt-near_cor) >  60:
#                                        cor_events[cor_count][0]=rt
#                                        cor_events[cor_count][1]="%1d" % p
#                                        cor_events[cor_count][2]="%1d" % yr
#                                        cor_events[cor_count][3]="%1d" % jd
#                                        cor_events[cor_count][4]="%1d" % hr
#                                        cor_events[cor_count][5]="%1d" % mi
#                                        cor_events[cor_count][6]="%1d" % se
#                                        cor_events[cor_count][7]="%1d" % ms
#                                        cor_events[cor_count][8]="%1d" % Events_today[ind][2]
#                                        cor_events[cor_count][9]="%1d" % Events_today[ind][3]
#                                        cor_events[cor_count][10]="%1d" % Events_today[ind][4]
#                                        cor_events[cor_count][11]="%1d" % Events_today[ind][5]
#                                        cor_events[cor_count][12]="%1d" % Events_today[ind][6]
#                                        cor_events[cor_count][13]="%1d" % Events_today[ind][7]
#                                        cor_events[cor_count][14]="%1d" % Events_today[ind][8]
#                                        cor_events[cor_count][15]="%1d" % Events_today[ind][9]
#                                        cor_events[cor_count][16]="%1d" % Events_today[ind][10]
#                                        cor_events[cor_count][17]="%1d" % Events_today[ind][11]
#                                        cor_events[cor_count][18]="%1d" % Events_today[ind][12]
#                                        cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        correlated_today.append(expt)
#                                        cor_count += 1
#                                        Number_of_events_today += 1
#                                    else:
#                                        cor_events[inde][18]="%1d" % 1
#
#                        else:
#                            cor_events[cor_count][0]= rt
#                            cor_events[cor_count][1]="%1d" % p
#                            cor_events[cor_count][2]="%1d" % yr
#                            cor_events[cor_count][3]="%1d" % jd
#                            cor_events[cor_count][4]="%1d" % hr
#                            cor_events[cor_count][5]="%1d" % mi
#                            cor_events[cor_count][6]="%1d" % se
#                            cor_events[cor_count][7]="%1d" % ms
#                            cor_events[cor_count][18]="%1d" % 1
#                            cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            correlated_today.append(expt)
#                            cor_count += 1
#                            Number_of_events_today += 1
#                        
#                        lastcorr=0
#                        add += window-step     +30    #skip 60s to avoid double catch
#
#            else:
#                add += window-step
#                lastcorr=0
#    else:
#        switch5s=0
#                
#    #%% Station 6 scan  
#    if trace6s.stats.npts > 201:
#        trace6s.detrend(type='linear')
#        trace6s.detrend(type='demean')
#        trace6s.filter("bandpass", freqmin=fmin,freqmax=fmax)
#    
#        av_amp6s=1100
#   
#        begin = trace6s.stats.starttime
#        window=int(len(trs_s6)/100)
#        add=0
#        lastcorr=0 
#        expt=0
#        dl=trace6s.stats.endtime - trace6s.stats.starttime    
#        for t in range(0,int(trace6s.stats.npts/100), step):
#            view_start = begin + t + add
#            view_end = view_start + window
#            if view_start > trace6s.stats.endtime - window :
#                break
#    
#            trc6s = trace6s.slice(starttime = view_start  , endtime= view_end)
#            if trc6s.stats.npts < 100:
#                break
#            trc6s.detrend(type='linear')
#            trc6s.detrend(type='demean')
#            # only look for things above noise level 
#            amp=abs(trc6s.max())
#            
#            if amp > av_amp6s:       
#                peak6,cf6,bwid6,bwid256 = freq_info(trc6s.data,trc6s.stats.starttime,trc6s.stats.endtime)      
#                if (1.25 < cf6 < 3) and (1.1 < peak6 < 3) and (0.1 < bwid6 < 2.5): 
#                         #correlate between data and earthquake envelopes 
#                    trc6s_e = obspy.signal.filter.envelope(trc6s.data) 
#                    top_v,top,corell = corel(trs_s6,trc6s_e,shift)
#                    
#                    if 0.0 in corell : # in the case of missing data - which will cause function to crash
#                        break            
# 
#                    expt_new = view_start - (top/100) +10
#
#                    if abs(top_v) > crite and top_v > lastcorr:
#                        
#                        lastcorr = top_v
#                        expt = expt_new                                    
#                        
#                    if lastcorr > top_v:
#                        file_ls06.write(str(expt))
#                        file_ls06.write("\n")
#                        
#                        rt=expt.timestamp	
#                        jd=expt.julday
#                        yr=expt.year
#                        mo=expt.month
#                        da=expt.day
#                        hr=expt.hour
#                        mi=expt.minute
#                        se=expt.second
#                        ms=expt.microsecond/10000
#                        near,ind=find_nearest(Events_today[:,0], rt)
#
#                        if num_active[p] > 1:
#                            if (lb_num_active[p] == 0) and (abs(rt-near) > 60):
#                                Events_today[find_count][0]=rt
#                                Events_today[find_count][1]=p
#                                Events_today[find_count][13]=1
#                                Events_today = np.lib.pad(Events_today, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                All_events.append(expt)
#                                find_count += 1
#                            if abs(rt-near) <= 20:
#                                Events_today[ind][13]=1
#                                if (lb_num_active[p] == 0) or ((lb_num_active[p] > 0) and ((Events_today[ind][2] == 1) or (Events_today[ind][4] == 1))):
#                                    near_cor,inde=find_nearest(cor_events[:,0], rt)
#                                    if abs(rt-near_cor) >  60:
#                                        cor_events[cor_count][0]=rt
#                                        cor_events[cor_count][1]="%1d" % p
#                                        cor_events[cor_count][2]="%1d" % yr
#                                        cor_events[cor_count][3]="%1d" % jd
#                                        cor_events[cor_count][4]="%1d" % hr
#                                        cor_events[cor_count][5]="%1d" % mi
#                                        cor_events[cor_count][6]="%1d" % se
#                                        cor_events[cor_count][7]="%1d" % ms
#                                        cor_events[cor_count][8]="%1d" % Events_today[ind][2]
#                                        cor_events[cor_count][9]="%1d" % Events_today[ind][3]
#                                        cor_events[cor_count][10]="%1d" % Events_today[ind][4]
#                                        cor_events[cor_count][11]="%1d" % Events_today[ind][5]
#                                        cor_events[cor_count][12]="%1d" % Events_today[ind][6]
#                                        cor_events[cor_count][13]="%1d" % Events_today[ind][7]
#                                        cor_events[cor_count][14]="%1d" % Events_today[ind][8]
#                                        cor_events[cor_count][15]="%1d" % Events_today[ind][9]
#                                        cor_events[cor_count][16]="%1d" % Events_today[ind][10]
#                                        cor_events[cor_count][17]="%1d" % Events_today[ind][11]
#                                        cor_events[cor_count][18]="%1d" % Events_today[ind][12]
#                                        cor_events[cor_count][19]="%1d" % Events_today[ind][13]
#                                        cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                                        correlated_today.append(expt)
#                                        cor_count += 1
#                                        Number_of_events_today += 1
#                                    else:
#                                        cor_events[inde][19]="%1d" % 1
#
#                        else:
#                            cor_events[cor_count][0]=rt
#                            cor_events[cor_count][1]="%1d" % p
#                            cor_events[cor_count][2]="%1d" % yr
#                            cor_events[cor_count][3]="%1d" % jd
#                            cor_events[cor_count][4]="%1d" % hr
#                            cor_events[cor_count][5]="%1d" % mi
#                            cor_events[cor_count][6]="%1d" % se
#                            cor_events[cor_count][7]="%1d" % ms
#                            cor_events[cor_count][19]="%1d" % 1
#                            cor_events = np.lib.pad(cor_events, ((0,1),(0,0)), 'constant', constant_values=(0))
#                            correlated_today.append(expt)
#                            cor_count += 1
#                            Number_of_events_today += 1
#                        lastcorr=0                                        
#                        add += window-step +30              #skip 60s to avoid double catch
#
#            else:
#                add += window-step
#                lastcorr=0 
#    else:
#        switch6s=0
#

               
#%% collect information from each station for single day and save in arrays                   
    Events_per_day[p][0]= Number_of_events_today
    Events_per_day[p][1]=p
    Events_per_day=np.lib.pad(Events_per_day, ((0,1),(0,0)), 'constant', constant_values=(0))    
    todays_event.sort()
    correlated_today.sort()
    for l in range(0,len(todays_event)):
        file_exp.write(str(todays_event[l]))
        file_exp.write("\n")
    
       
    for s in range(0,len(correlated_today)):
        file_cor2.write(str(correlated_today[s]))
        file_cor2.write("\n")       

#%% 
print("End of Scan")
np.savetxt("/Users/william/Documents/scanner/all_stations/Explosions_per_day_LB03_2018.csv", Events_per_day,delimiter=",",header="Number_of_Events,Day_number")            

col=0       
cor_events=cor_events[np.argsort(cor_events[:,col])] 
if cor_events[0,0]== 0.0:
    cor_events = np.delete(cor_events, (0), axis=0)       
np.savetxt("/Users/william/Documents/scanner/all_stations/EXP_all_coincidence_LB03_2018.csv", cor_events,delimiter=",",header="Time_stamp,Day_number,Year,Day_of_year,Hour,Min,Sec,Milisec,LBO1,LB02,LB03,LB04,LB05,LB06,LSO1,LS02,LS03,LS04,LS05,LS06")

data_stream = Stream()
for r in range(0,len(cor_events)):
#    if cor_events[r][8]==1:
#        tr1a=stream1[int(cor_events[r][1])]
#        dstart=UTCDateTime(cor_events[r][0])
#        tr1_save = tr1a.slice(starttime=dstart-15, endtime=dstart+65)
#        tr1_save.detrend(type='demean')
#        tr1_save.detrend(type='linear')
#        if type(tr1_save[0]) == np.float64:
#            data_stream.append(tr1_save)         
#    if cor_events[r][9]==1:
#        tr2a=stream2[int(cor_events[r][1])]
#        dstart=UTCDateTime(cor_events[r][0])
#        tr2_save = tr2a.slice(starttime=dstart-15, endtime=dstart+65)
#        tr2_save.detrend(type='demean')
#        tr2_save.detrend(type='linear')
#        if type(tr2_save[0]) == np.float64:
#            data_stream.append(tr2_save)
    if cor_events[r][10]==1:
        tr3a=stream3[int(cor_events[r][1])]
        dstart=UTCDateTime(cor_events[r][0])
        tr3_save = tr3a.slice(starttime=dstart-15, endtime=dstart+65)
        tr3_save.detrend(type='demean')
        tr3_save.detrend(type='linear')
        if type(tr3_save[0]) == np.float64:
            data_stream.append(tr3_save)         
#    if cor_events[r][11]==1:
#        tr4a=stream4[int(cor_events[r][1])]
#        dstart=UTCDateTime(cor_events[r][0])
#        tr4_save = tr4a.slice(starttime=dstart-15, endtime=dstart+65)
#        tr4_save.detrend(type='demean')
#        tr4_save.detrend(type='linear')
#        if type(tr4_save[0]) == np.float64:
#            data_stream.append(tr4_save)
#    if cor_events[r][12]==1:
#        tr5a=stream5[int(cor_events[r][1])]
#        dstart=UTCDateTime(cor_events[r][0])
#        tr5_save = tr5a.slice(starttime=dstart-15, endtime=dstart+65)
#        tr5_save.detrend(type='demean')
#        tr5_save.detrend(type='linear')
#        if type(tr5_save[0]) == np.float64:
#            data_stream.append(tr5_save)
#    if cor_events[r][13]==1:
#        tr6a=stream6[int(cor_events[r][1])] 
#        dstart=UTCDateTime(cor_events[r][0])
#        tr6_save = tr6a.slice(starttime=dstart-15, endtime=dstart+65)
#        tr6_save.detrend(type='demean')
#        tr6_save.detrend(type='linear')
#        if type(tr6_save[0]) == np.float64:
#            data_stream.append(tr6_save)
#    if cor_events[r][14]==1:
#        tr1sa=stream1s[int(cor_events[r][1])]
#        dstart=UTCDateTime(cor_events[r][0])
#        tr1s_save = tr1sa.slice(starttime=dstart-15, endtime=dstart+65)
#        tr1s_save.detrend(type='demean')
#        tr1s_save.detrend(type='linear')
#        if type(tr1s_save[0]) == np.float64:
#            data_stream.append(tr1s_save)         
#    if cor_events[r][15]==1:
#        tr2sa=stream2s[int(cor_events[r][1])]
#        dstart=UTCDateTime(cor_events[r][0])
#        tr2s_save = tr2sa.slice(starttime=dstart-15, endtime=dstart+65)
#        tr2s_save.detrend(type='demean')
#        tr2s_save.detrend(type='linear')
#        if type(tr2s_save[0]) == np.float64:
#            data_stream.append(tr2s_save)
#    if cor_events[r][16]==1:
#        tr3sa=stream3s[int(cor_events[r][1])]
#        dstart=UTCDateTime(cor_events[r][0])
#        tr3s_save = tr3sa.slice(starttime=dstart-15, endtime=dstart+65)
#        tr3s_save.detrend(type='demean')
#        tr3s_save.detrend(type='linear')
#        if type(tr3s_save[0]) == np.float64:
#            data_stream.append(tr3s_save)         
#    if cor_events[r][17]==1:
#        tr4sa=stream4s[int(cor_events[r][1])]
#        dstart=UTCDateTime(cor_events[r][0])
#        tr4s_save = tr4sa.slice(starttime=dstart-15, endtime=dstart+65)
#        tr4s_save.detrend(type='demean')
#        tr4s_save.detrend(type='linear')
#        if type(tr4s_save[0]) == np.float64:
#            data_stream.append(tr4s_save)
#    if cor_events[r][18]==1:
#        tr5sa=stream5s[int(cor_events[r][1])]
#        dstart=UTCDateTime(cor_events[r][0])
#        tr5s_save = tr5sa.slice(starttime=dstart-15, endtime=dstart+65)
#        tr5s_save.detrend(type='demean')
#        tr5s_save.detrend(type='linear')
#        if type(tr5s_save[0]) == np.float64:
#            data_stream.append(tr5s_save)
#    if cor_events[r][19]==1:
#        tr6sa=stream6s[int(cor_events[r][1])] 
#        dstart=UTCDateTime(cor_events[r][0])
#        tr6s_save = tr6sa.slice(starttime=dstart-15, endtime=dstart+65)
#        tr6s_save.detrend(type='demean')
#        tr6s_save.detrend(type='linear')
#        if type(tr6s_save[0]) == np.float64:
#            data_stream.append(tr6s_save)

if len(data_stream) > 0:            
    data_stream.write('/Users/william/Documents/scanner/output_data/EXP_LB03_2018', format='MSEED')
#%%

file_exp.close() 
file_cor2.close()
#file_lb01.close()  
#file_lb02.close()  
file_lb03.close()  
#file_lb04.close()  
#file_lb05.close()  
#file_lb06.close()  
#file_ls01.close()  
#file_ls02.close()  
#file_ls03.close()  
#file_ls04.close()  
#file_ls05.close()  
#file_ls06.close()           
print("Data saved") 
#T2=time.clock()
#elapsed= T2-T1  
#print('Time taken:', elapsed)   
#    
    
    
    
    
    