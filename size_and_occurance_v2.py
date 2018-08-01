#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 18:39:00 2018

@author: william
"""

import obspy
from obspy import read
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from obspy import Stream
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import UTCDateTime


#%% import data
stre = read("/Users/william/Documents/scanner/output_data/EXP_all_data_stream_month_1.mseed")
per_day = genfromtxt("/Users/william/Documents/scanner/all_stations/Explosions_per_day_v2.csv", delimiter=',',skip_header=1,skip_footer=1)
data = genfromtxt("/Users/william/Documents/scanner/all_stations/EXP_all_coincidence_month_1.csv", delimiter=',',skip_header=1)

#%% Events per day and per week
plt.figure(10001)
plt.plot(per_day[:,1],per_day[:,0])
plt.ylim([0,max(per_day[:,0])+10])
plt.xlabel('Day Number')
plt.ylabel('Number of Explosions')
plt.title('Explosions per Day')

epw=np.zeros(shape=(1,2))
day=0
week=0
nepw=0
for x in range(0,len(per_day)):
    nepw += per_day[x,0]
    day += 1
    if day == 7:
        epw[week][0] = nepw
        epw[week][1] = week
        week +=1
        day=0
        nepw=0
        if len(per_day)-x > 7:
            epw = np.lib.pad(epw, ((0,1),(0,0)), 'constant', constant_values=(0))
plt.figure(10002)
plt.plot(epw[:,1],epw[:,0])
plt.ylim([0,max(epw[:,0])+50])
plt.xlabel('Week Number')
plt.ylabel('Number of Explosions')
plt.title('Explosions per Week')
#

#%%    Seismic calibrations

#Broadband
#LB01 is redeployed with different calibration factors
LB01sc1=0.000001/750            # before 2015-12-05T00:00:01.000000Z
LB01sc2=(10*0.000001)/(750*256) # 2015-12-05T00:00:01 - 2016-06-15T00:00:01 + maybe till end

LB02sc=0.000001/750
LB03sc=0.000001/(750*256)
LB04sc=0.000001/(750*256)
LB05sc=0.000001/(750*256)
LB06sc=(10*0.000001)/750

# Short period
LSsc=0.000000122/800

#%%    Acoustic calibrations

LB01ac = 0.000001/0.0250
LB02ac = 0.000001/0.0250
LB03ac = 0.000001/(0.0250*256)
LB04ac = 0.000001/(0.0250*256)
LB05ac = 0.000001/(0.0250*256)

#%%   Station distances



#%% constants

rho_atmos = 1.2 # kg/m3 at 15˚C
c_atmos = 340 # m/s at 15˚C

rho_earth = 2000 # kg/m3
c_earth = 2500 # m/s 

#%%





#%% "Energy" of each event

#################### Need to use proper equations with calibrations for all below #######################

event_stream = Stream()
event_list=np.zeros(shape=(1,1))
event_count=0

event_stream.append(stre[0])
event_list[event_count]=stre[0].stats.starttime.timestamp
event_list = np.lib.pad(event_list, ((0,1),(0,0)), 'constant', constant_values=(0))
event_count +=1                            
for x in range(1,len(stre)):
    if stre[x].stats.station == "LB01":
        event_stream.append(stre[x])
        event_list[event_count]=stre[x].stats.starttime.timestamp
        event_list = np.lib.pad(event_list, ((0,1),(0,0)), 'constant', constant_values=(0))
        event_count +=1 
    else:
        rt=stre[x].stats.starttime.timestamp
        near,ix=find_nearest(event_list[:,0], rt)
        if abs(near-rt) > 60:
            event_stream.append(stre[x])
            event_list[event_count]=stre[x].stats.starttime.timestamp
            event_list = np.lib.pad(event_list, ((0,1),(0,0)), 'constant', constant_values=(0))
            event_count +=1 


#print(len(data_stream))

sr = 100
nsta=int(1*sr)                                      
nlta=int(10*sr)                                     
trig_on=2.5                                          
trig_off=0.05

#for x in range(0,len(data_stream)):
for x in range(0,10):
    data_s=event_stream[x].data
    max_a = data_s.max()
    min_a = data_s.min()
    p2p= max_a-min_a
    cft=recursive_sta_lta(data_s, nsta, nlta)
#    plot_trigger(sq_stream[x], cft, trig_on, trig_off)     
    on_off = trigger_onset(cft,trig_on,trig_off)                                     
    start = event_stream[x].stats.starttime
    tr = event_stream[x].slice(starttime=start+(on_off[0,0]/sr) , endtime=start+(on_off[0,1]/sr)) 
    print('event:',event_stream[x].stats.starttime ,'from station: ',stre[x].stats.station,', has energy: ',sum(np.square(tr.data)),' and peak to peak:', p2p)
    plt.figure(x)
    plt.plot(tr) 
    plt.figure(x+20)
    plt.plot(np.square(tr.data))   
     
#%% energy information
    
        
day_one = 1416787200.0
day_energy=0
day_count=0
e_count=0
days=0
energy_each_event = av_energy_list =np.zeros(shape=(1,1))
av_energy_list =np.zeros(shape=(1,1))
total_energy_list =np.zeros(shape=(1,1))

for x in range(0,len(per_day)):#len(per_day)
    day_start = day_one + x*24*60*60
    day_end = day_start + 24*60*60 - 0.01
#    print(UTCDateTime(day_start),' to', UTCDateTime(day_end))
    for p in range(0,len(event_stream)):
        if day_start < event_stream[p].stats.starttime.timestamp < day_end :
            # ADD IN CALIBRATIONS #
            day_count += 1
            data_s=event_stream[p].data
            max_a = data_s.max()
            min_a = data_s.min()
            p2p= max_a-min_a
            cft=recursive_sta_lta(data_s, nsta, nlta)
        #    plot_trigger(sq_stream[x], cft, trig_on, trig_off)     
            on_off = trigger_onset(cft,trig_on,trig_off)                                     
            start = event_stream[p].stats.starttime
            tr = event_stream[p].slice(starttime=start+(on_off[0,0]/sr) , endtime=start+(on_off[0,1]/sr)) 
            event_energy = sum(np.square(tr.data))
            day_energy += event_energy
#            print('event:',event_stream[p].stats.starttime ,'station: ',event_stream[p].stats.station,', energy: ',event_energy,', peak to peak:', p2p)
            energy_each_event[e_count] = event_energy
            energy_each_event = np.lib.pad(energy_each_event, ((0,1),(0,0)), 'constant', constant_values=(0))
            e_count += 1
            
#            if event_energy > 1e16:
#                tr.plot()
            
    print('total energy in the day:',UTCDateTime(day_start),'=', day_energy)
    av_day_energy = day_energy/day_count
    print('average energy in day:',UTCDateTime(day_start),'=', av_day_energy)
    
    av_energy_list[days]= av_day_energy
    total_energy_list[days] = day_energy 
    av_energy_list = np.lib.pad(av_energy_list, ((0,1),(0,0)), 'constant', constant_values=(0))
    total_energy_list = np.lib.pad(total_energy_list, ((0,1),(0,0)), 'constant', constant_values=(0))
    day_count = 0
    day_energy=0
    days += 1
    
    


plt.figure(2001)
plt.plot(av_energy_list)

plt.figure(2002)
plt.plot(total_energy_list)

plt.figure(2003)
plt.plot(energy_each_event)























