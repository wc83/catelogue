#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 09:33:56 2018

@author: william
"""

def stacked_LBz_exp():    
    import numpy as np
    import obspy
    from obspy.core import read   
    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    from obspy import Stream
    import matplotlib.pyplot as plt
    #%%
    
    seisb = read('/Users/william/Documents/lb01/14_326z.mseed')
    stream1=seisb.copy()
    stream1=stream1.clear()
    stream2=stream1.copy()
    stream3=stream1.copy()
    stream4=stream1.copy()
    stream5=stream1.copy()
    stream6=stream1.copy()
    stream7=stream1.copy()
    stream1b=stream1.copy()
    stream2b=stream1.copy()
    stream3b=stream1.copy()
    stream4b=stream1.copy()
    stream5b=stream1.copy()
    stream6b=stream1.copy()
    stream7b=stream1.copy()
    
    year1=2014
    month1=11
    day1=22
    hour1=0
    minute1=0
    second1=0
    fmin=0.1
    fmax=10
    #%% LB01
    sta = 'LB01' # STATION LB01
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(0,10,1):
        if x==0:
            M=[[1,36,46],[2,37,17],[4,44,34],[12,19,2],[20,21,27]]
        if x==1:
            M=[[3,3,40],[5,59,54],[9,47,32],[15,18,34],[18,39,4]]
        if x==2:
            M=[[2,50,57],[8,22,15],[13,11,9],[19,7,58],[22,52,31]]
        if x==3:
            M=[[1,49,8],[2,25,56],[7,38,37],[8,56,5],[13,40,50]]
        if x==4:
            M=[[2,22,9],[3,14,20],[7,6,43],[11,6,13],[17,36,1]]
        if x==5:
            M=[[1,25,18],[12,59,20],[15,5,46],[19,11,43],[21,33,7]]
        if x==6:
            M=[[3,21,27],[4,20,18],[7,7,38],[12,42,54],[22,16,3]]
        if x==7:
            M=[[3,23,12],[4,6,44],[5,49,59],[15,51,20],[19,46,30]]
        if x==8:
            M=[[1,9,2],[5,30,42],[12,31,5],[20,15,40],[22,27,20]]
        if x==9:
            M=[[1,22,40],[5,27,21],[10,13,51],[12,41,20],[22,1,41]]
            
        for i in range(0,len(M),1):   
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]
            
            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  - 40
            t2 = t1 + 100 
            seis = Stream()
            seis = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc1 = seis[0].slice(starttime = t1 + 30  , endtime= t2 - 10)
            trc1.detrend(type='demean')
            trc1.detrend(type='linear')
            stream1.append(trc1)    
    #        trc.plot(type='relative',color='b')
    
    for x in range(0,len(stream1)):
        trc=stream1[x].normalize()
        stream1b.append(trc)
    
    stack_norm1 = np.sum([abs(trc.data) for trc in stream1b], axis=0)
    stack_norm1 = stack_norm1/len(stream1b)
    
#    plt.figure(1)
#    plt.plot(stack_norm1,color='r')
#    plt.title('LB01 stacked explosion waveform')    
#        
    #%% LB02
    sta = 'LB02' # STATION LB02
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(0,10,1):
        if x==0:
            M=[[1,36,46],[2,37,17]]
        if x==1:
            M=[[5,59,54],[9,47,32],[18,39,4]]
        if x==2:
            M=[[13,11,9],[22,52,31]]
        if x==3:
            M=[[1,49,8],[2,25,56],[8,56,5]]
        if x==4:
            M=[[2,22,9],[7,6,43],[11,6,13]]
        if x==5:
            M=[[1,25,18],[12,59,20],[15,5,46],[19,11,42]]
        if x==6:
            M=[[4,20,18],[7,7,38]]
        if x==7:
            M=[[3,23,12],[5,49,59],[15,51,20]]
        if x==8:
            M=[[12,31,5],[22,27,20]]
        if x==9:
            M=[[1,22,40],[12,41,23]]
            
        for i in range(0,len(M),1):   
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]
            
            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -40
            t2 = t1 + 100
            seis = Stream()
            seis = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc2 = seis[0].slice(starttime = t1 + 30  , endtime= t2 - 10)
            trc2.detrend(type='demean')
            trc2.detrend(type='linear')
            stream2.append(trc2)    
    #        trc.plot(type='relative',color='b')
    
    for x in range(0,len(stream2)):
        trc=stream2[x].normalize()
        stream2b.append(trc)
    
    stack_norm2 = np.sum([abs(trc.data) for trc in stream2b], axis=0)
    stack_norm2 = stack_norm2/len(stream2b)
#    
#    plt.figure(2)
#    plt.plot(stack_norm2,color='r')
#    plt.title('LB02 stacked explosion waveform')    
            
    #%% LB03
    sta = 'LB03' # STATION LB03
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(2,10,1):
        if x==2:
            M=[[2,50,57],[8,22,15],[13,11,8],[19,7,58],[22,52,31]]
        if x==3:
            M=[[1,49,8],[2,25,56],[7,38,37],[8,56,5],[13,40,50]]
        if x==4:
            M=[[2,22,9],[3,14,20],[7,6,43],[11,6,13],[17,36,1]]
        if x==5:
            M=[[1,25,18],[12,59,20],[15,5,46],[19,11,42],[21,33,7]]
        if x==6:
            M=[[3,21,27],[12,42,54],[22,16,0]]
        if x==7:
            M=[[3,23,14],[4,6,44],[15,51,19],[19,46,30]]
        if x==8:
            M=[[1,9,2],[5,30,42],[12,31,5],[20,15,40],[22,27,19]]
        if x==9:
            M=[[1,22,40],[5,27,21],[10,13,51],[12,41,23],[22,1,41]]
        for i in range(0,len(M),1):   
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]
            
            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -40
            t2 = t1 + 100
            seis = Stream()
            seis = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc3 = seis[0].slice(starttime = t1 + 30  , endtime= t2 - 10)
            trc3.detrend(type='demean')
            trc3.detrend(type='linear')
            stream3.append(trc3)    
    #        trc3.plot(type='relative',color='b')
    
    for x in range(0,len(stream3)):
        trc=stream3[x].normalize()
        stream3b.append(trc)
    
    stack_norm3 = np.sum([abs(trc.data) for trc in stream3b], axis=0)
    stack_norm3 = stack_norm3/len(stream3b)
    
#    plt.figure(3)
#    plt.plot(stack_norm3,color='r')
#    plt.title('LB03 stacked explosion waveform')      
        
    #%% LB04
    sta = 'LB04' # STATION LB04
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(3,10,1):
        if x==3:
            M=[[1,49,8],[2,25,57],[7,38,37],[8,56,5],[13,40,51]]
        if x==4:
            M=[[2,22,11],[3,14,21],[7,6,43],[11,6,16],[17,36,1]]
        if x==5:
            M=[[1,25,20],[12,59,20],[15,5,46],[19,11,43],[21,33,9]]
        if x==6:
            M=[[3,21,27],[4,20,18],[7,7,38],[12,42,54],[22,16,1]]
        if x==7:
            M=[[3,23,14],[4,6,44],[5,50,0],[15,51,19],[19,46,30]]
        if x==8:
            M=[[1,9,2],[5,30,42],[12,31,5],[20,15,40],[22,27,20]]
        if x==9:
            M=[[1,22,41],[5,27,21],[10,13,51],[12,41,25],[22,1,43]]
            
        for i in range(0,len(M),1):   
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]
            
            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -40
            t2 = t1 +100
            seis = Stream()
            seis = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc4 = seis[0].slice(starttime = t1 + 30  , endtime= t2 - 10)
            trc4.detrend(type='demean')
            trc4.detrend(type='linear')
            stream4.append(trc4)    
    #        trc.plot(type='relative',color='b')
    
    for x in range(0,len(stream4)):
        trc=stream4[x].normalize()
        stream4b.append(trc)
    
    stack_norm4 = np.sum([abs(trc.data) for trc in stream4b], axis=0)
    stack_norm4 = stack_norm4/len(stream4b)
    
#    plt.figure(4)
#    plt.plot(stack_norm4,color='r')
#    plt.title('LB04 stacked explosion waveform')      
        
    #%% LB05
    sta = 'LB05' # STATION LB05
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(4,10,1):
        if x==4:
            M=[[2,22,11],[3,14,22],[7,6,44],[11,6,14],[17,36,0]]
        if x==5:
            M=[[1,25,18],[12,59,20],[15,5,46],[19,11,42],[21,33,8]]
        if x==6:
            M=[[7,7,38],[12,42,54],[22,16,1]]
        if x==7:
            M=[[3,23,14],[4,6,44],[5,50,1]]
        if x==8:
            M=[[1,9,3],[5,30,42],[12,31,5],[20,15,42],[22,27,20]]
        if x==9:
            M=[[1,22,42],[5,27,21],[10,13,51],[12,41,24],[22,1,42]]
            
        for i in range(0,len(M),1):   
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]
            
            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -40
            t2 = t1 +100
            seis5 = Stream()
            seis5 = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis5[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc5 = seis5[0].slice(starttime = t1 + 30  , endtime= t2 - 10)
            trc5.detrend(type='demean')
            trc5.detrend(type='linear')
            stream5.append(trc5)    
    #        trc.plot(type='relative',color='b') 

    for x in range(0,len(stream5)):
        trc=stream5[x].normalize()
        stream5b.append(trc)
    
    stack_norm5 = np.sum([abs(trc.data) for trc in stream5b], axis=0)
    stack_norm5 = stack_norm5/len(stream5b)
    
#    plt.figure(5)
#    plt.plot(stack_norm5,color='r')
#    plt.title('LB05 stacked explosion waveform')  
#    
    #%% LB06    
     
    # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
    sta = 'LB06' # STATION LB01
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(2016, 6, 21, hour1, minute1, second1) #the format is year:day_of_the_year:month
     
    M=[[4,40,32,0],[16,2,38,8],[4,42,4,22],[5,14,34,22],[6,43,15,33],[10,59,20,33],[1,11,17,36]]
    
    for i in range(0,len(M)):
        h=M[i][0]
        m=M[i][1]
        s=M[i][2]
        d=M[i][3]
        
        t1 = t0 + d*24*60*60 + h*60*60 + m*60 + s  -60
        t2 = t1 +110
        seis6 = Stream()
        seis6 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        seis6[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
        trc6 = seis6[0].slice(starttime = t1 + 30  , endtime= t2 - 10)
        trc6.detrend(type='demean')
        trc6.detrend(type='linear')
        stream6.append(trc6)    
    #    trc6.plot(type='relative',color='b')
    
    for x in range(0,len(stream6)):
        trc=stream6[x].normalize()
        stream6b.append(trc)
    
    stack_norm6 = np.sum([abs(trc.data) for trc in stream6b], axis=0)
    stack_norm6 = stack_norm6/len(stream6)
    
#    plt.figure(6)
#    plt.plot(stack_norm6,color='r')
#    plt.title('LB06 stacked explosion waveform') 
        
    #stack_norm6b = obspy.signal.filter.envelope(trc6.data)
    
    #%% LB07   
#     
#    # STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
#    sta = 'LB07' # STATION LB01
#    cha = 'HHZ' # CHANNEL - Vertical
#    net = 'Z4'  # Santiaguito volcano
#    loc = ''    # location, it depends mostly of which network you are in. 
#    
#    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#    t0 = UTCDateTime(2016, 6, 21, hour1, minute1, second1) #the format is year:day_of_the_year:month
#     
#    M=[[5,14,33,22],[10,59,27,33],[1,11,29,36]]
#    
#    for i in range(0,len(M)):
#        h=M[i][0]
#        m=M[i][1]
#        s=M[i][2]
#        d=M[i][3]
#        
#        t1 = t0 + d*24*60*60 + h*60*60 + m*60 + s  -60
#        t2 = t1 +110
#        seis7 = Stream()
#        seis7 = client.get_waveforms(net, sta, '', cha, t1 , t2)
#        
#        seis7[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
#        trc7 = seis7[0].slice(starttime = t1 + 30  , endtime= t2 - 10)
#        trc7.detrend(type='demean')
#        trc7.detrend(type='linear')
#        stream7.append(trc7)    
##        trc7.plot(type='relative',color='b')
#    
#    for x in range(0,len(stream7)):
#        trc=stream7[x].normalize()
#        stream7b.append(trc)
#    
#    stack_norm7 = np.sum([abs(trc.data) for trc in stream7b], axis=0)
#    stack_norm7 = stack_norm7/len(stream7)
#    
##    plt.figure(7)
##    plt.plot(stack_norm7,color='r') 
##    plt.title('LB07 stacked explosion waveform')  
#    
#    #stack_norm7b = obspy.signal.filter.envelope(trc7.data)
    #%%
    return(stack_norm1,stack_norm2,stack_norm3,stack_norm4,stack_norm5,stack_norm6)#,stack_norm7)
