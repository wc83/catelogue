#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 16:44:01 2018

@author: william
"""


def get_LBz(day):
    from obspy.core import read
        
    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    
    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    #from scipy.signal import welch
    from obspy import Stream
    
    year1=2014
    month1=11
    day1=24
    hour1=0
    minute1=0
    second1=0
#%% LB01    
    try:  
        
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB01' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)

    except: # give 2 seconds of data instead
        
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB01' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 
        t2 = t1 + 2 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, '', cha, t1 , t2)

        st1.detrend(type='linear')
        st1.detrend(type='demean')
        st1[0].filter("bandpass", freqmin=0.1,freqmax=0.1)
    
#%%    LB02

    try:  
     
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB02' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
        
        client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st2 = Stream()
        st2 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st2.detrend(type='linear')
        st2.detrend(type='demean')
        break_test=st2
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)

    except: # give 2 seconds of data instead

        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB01' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 
        t2 = t1 + 2 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st2 = Stream()
        st2 = client.get_waveforms(net, sta, '', cha, t1 , t2)

        st2.detrend(type='linear')
        st2.detrend(type='demean')
        st2[0].filter("bandpass", freqmin=0.1,freqmax=0.1)
  
    
    
#%% LB03    
    try:  
     
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB03' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st3 = Stream()
        st3 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st3.detrend(type='linear')
        st3.detrend(type='demean')
        break_test=st3
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)

    except: # give 2 seconds of data instead
        
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB01' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 
        t2 = t1 + 2 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st3 = Stream()
        st3 = client.get_waveforms(net, sta, '', cha, t1 , t2)

        st3.detrend(type='linear')
        st3.detrend(type='demean')
        st3[0].filter("bandpass", freqmin=0.1,freqmax=0.1)
  
#%% LB04    
    try:  
     
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB04' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st4 = Stream()
        st4 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st4.detrend(type='linear')
        st4.detrend(type='demean')
        break_test=st4
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)

    except: # give 2 seconds of data instead
        
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB01' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 
        t2 = t1 + 2 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st4 = Stream()
        st4 = client.get_waveforms(net, sta, '', cha, t1 , t2)

        st4.detrend(type='linear')
        st4.detrend(type='demean')
        st4[0].filter("bandpass", freqmin=0.1,freqmax=0.1)    
    
#%% LB05    
    try:  
     
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB05' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st5 = Stream()
        st5 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st5.detrend(type='linear')
        st5.detrend(type='demean')
        break_test=st5
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)

    except: # give 2 seconds of data instead
        
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB01' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 
        t2 = t1 + 2 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st5 = Stream()
        st5 = client.get_waveforms(net, sta, '', cha, t1 , t2)

        st5.detrend(type='linear')
        st5.detrend(type='demean')
        st5[0].filter("bandpass", freqmin=0.1,freqmax=0.1)
    
#%% LB06    
    try:  
     
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB06' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st6 = Stream()
        st6 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st6.detrend(type='linear')
        st6.detrend(type='demean')
        break_test=st6
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)

    except: # give 2 seconds of blank data instead
        
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB01' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 
        t2 = t1 + 2 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st6 = Stream()
        st6 = client.get_waveforms(net, sta, '', cha, t1 , t2)

        st6.detrend(type='linear')
        st6.detrend(type='demean')
        st6[0].filter("bandpass", freqmin=0.1,freqmax=0.1)
    
    
#%% LB07
    try:  
     
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB07' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 
               
        client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st7 = Stream()
        st7 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st7.detrend(type='linear')
        st7.detrend(type='demean')
        break_test=st7
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)

    except: # give 2 seconds of blank data instead
        
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB01' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 
        t2 = t1 + 2 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st7 = Stream()
        st7 = client.get_waveforms(net, sta, '', cha, t1 , t2)

        st7.detrend(type='linear')
        st7.detrend(type='demean')
        st7[0].filter("bandpass", freqmin=0.1,freqmax=0.1)
    #%% return all stations
    return(st1,st2,st3,st4,st5,st6,st7)