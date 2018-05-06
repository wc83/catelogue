#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 10:24:41 2018

@author: william
"""


def get_LB01z(day):
    year1=2014
    month1=11
    day1=24
    hour1=0
    minute1=0
    second1=0
    
    from obspy.core import read
        
    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    
    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    #from scipy.signal import welch
    from obspy import Stream
    
    # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
    sta = 'LB01' # STATION LB01
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    # Corner frequency for high-pass filter
    #    hp_corner = 0.05
    
    # t1. and t2 are in hours:minutes:seconds
    # Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist 
    try:     
        client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 + day*24*60*60
        t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        #st.plot(method="full")
        st.detrend(type='linear')
        st.detrend(type='demean')
        break_test=st
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        #st.plot(color='b',starttime=t1, endtime=t2)
    #    print(st)
    except IndexError:
        year1=2014
        month1=11
        day1=24
        hour1=0
        minute1=0
        second1=0
        
        from obspy.core import read
            
        from obspy.clients.earthworm import Client
        from obspy import UTCDateTime
        
        from obspy.clients.earthworm import Client
        from obspy import UTCDateTime
        #from scipy.signal import welch
        from obspy import Stream
        
        # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
        sta = 'LB03' # STATION LB01
        cha = 'HHZ' # CHANNEL - Vertical
        net = 'Z4'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
        t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
        t1 = t0 
        t2 = t1 + 2 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        st = Stream()
        st = client.get_waveforms(net, sta, '', cha, t1 , t2)

        st.detrend(type='linear')
        st.detrend(type='demean')
        st[0].filter("bandpass", freqmin=0.1,freqmax=0.1)
    
    return(st)
    
    
    
    
    
    
    
    
    
    
    
    