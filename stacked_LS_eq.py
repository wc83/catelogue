#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 14:06:15 2018

@author: william
"""


def stacked_LSz_eq(fmi,fma):
    import numpy as np
    import obspy
    from obspy.core import read   
    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    from obspy import Stream
    import matplotlib.pyplot as plt
    #%%
    
    seisb = Stream()
    stream1=seisb.copy()
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
    day1=24
    hour1=0
    minute1=0
    second1=0
    fmin=fmi
    fmax=fma
    #%% LB01
    sta = 'LS01' # STATION LS01
    cha = 'EHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(0,13):
        if x==0:
            M=[[1,57,21],[18,42,57.5]]
        if x==1:
            M=[[9,16,54],[13,12,11]]
        if x==2:
            M=[[3,52,31],[5,14,25],[8,13,11]]
        if x==3:
            M=[[8,4,17.5]]
        if x==4:
            M=[]
        if x==5:
            M=[[7,56,47]]
        if x==6:
            M=[[12,47,4.5]]
        if x==7:
            M=[[8,37,51]]
        if x==8:
            M=[[9,45,59]]
        if x==9:
            M=[[10,34,43]]  
        if x==10:
            M=[[6,5,26.5]]  
        if x==11:
            M=[]  
        if x==12:
            M=[[14,25,21]]  
    #       
        for i in range(0,len(M),1):   
            
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]
    
            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -20
            t2 = t1 + 80
            seis1 = Stream()
            seis1 = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis1[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc1 = seis1[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
            trc1.detrend(type='demean')
            trc1.detrend(type='linear')
            stream1.append(trc1)
            trc1.plot(type='relative',color='b')
    
    
    for x in range(0,len(stream1)):
        trc=stream1[x].normalize()
        stream1b.append(trc)
    
    stack_norm1 = np.sum([abs(trc.data) for trc in stream1b], axis=0)
    stack_norm1 = stack_norm1/len(stream1b)
    
    plt.figure(1)
    plt.plot(stack_norm1,color='r')
    plt.title('LS01 stacked earthquake waveform')
    
    #%% LB02
    sta = 'LS02' # STATION LB02
    cha = 'EHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(0,13,1):
        if x==0:
            M=[[1,57,20],[8,30,16],[18,42,56.5],[18,59,58.5]]
        if x==1:
            M=[]
        if x==2:
            M=[[3,52,30],[5,14,24],[6,1,34],[8,13,9]]
        if x==3:
            M=[[6,50,48]]
        if x==4:
            M=[]
        if x==5:
            M=[[7,56,46]]
        if x==6:
            M=[[12,47,5],[12,55,3]]
        if x==7:
            M=[[8,37,51],[12,17,50.5]]
        if x==8:
            M=[[9,46,0]]
        if x==9:
            M=[[10,34,43.5]]  
        if x==10:
            M=[[14,53,22],[20,3,3.5]]  
        if x==11:
            M=[[19,58,46.5]]  
        if x==12:
            M=[[8,59,50],[14,25,21]]  
           
        for i in range(0,len(M),1):   
            
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]
    
            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -20
            t2 = t1 + 80
            seis2 = Stream()
            seis2 = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis2[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc2 = seis2[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
            trc2.detrend(type='demean')
            trc2.detrend(type='linear')
            stream2.append(trc2)
            trc2.plot(type='relative',color='b')
    
    for x in range(0,len(stream2)):
        trc=stream2[x].normalize()
        stream2b.append(trc)
    
    stack_norm2 = np.sum([abs(trc.data) for trc in stream2b], axis=0)
    stack_norm2 = stack_norm2/len(stream2b)
    
    plt.figure(2)
    plt.plot(stack_norm2,color='r')
    plt.title('LS02 stacked earthquake waveform')
    
    #%% LB03
    sta = 'LS03' # STATION LB02
    cha = 'EHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(1,13,1):
    
        if x==1:
            M=[[13,12,9.5]]
        if x==2:
            M=[[3,52,29.5],[5,14,23.5],[8,13,11]]
        if x==3:
            M=[]
        if x==4:
            M=[]
        if x==5:
            M=[[7,56,47]]
        if x==6:
            M=[[12,47,4],[12,55,2]]
        if x==7:
            M=[[8,37,49],[18,24,25]]
        if x==8:
            M=[]
        if x==9:
            M=[[10,34,41]]  
        if x==10:
            M=[[14,53,23.5],[20,3,5]]  
        if x==11:
            M=[[14,30,54],[19,58,48]]  
        if x==12:
            M=[[14,25,19.5]]  
           
        for i in range(0,len(M),1):   
            
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]
    
            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -20 
            t2 = t1 + 80  
            seis3 = Stream()
            seis3 = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis3[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc3 = seis3[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
            trc3.detrend(type='demean')
            trc3.detrend(type='linear')
            stream3.append(trc3)
            trc3.plot(type='relative',color='b')
    
    
    for x in range(0,len(stream3)):
        trc=stream3[x].normalize()
        stream3b.append(trc)
    
    stack_norm3 = np.sum([abs(trc.data) for trc in stream3b], axis=0)
    stack_norm3 = stack_norm3/len(stream3b)
    
    plt.figure(3)
    plt.plot(stack_norm3,color='r')
    plt.title('LS03 stacked earthquake waveform')
    
    #%% LB04
    sta = 'LS04' # STATION LB02
    cha = 'EHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(1,13,1):
    
        if x==1:
            M=[[9,16,53],[11,24,11],[13,12,8]]
        if x==2:
            M=[[3,52,28],[5,14,22],[6,1,31],[8,13,9],[9,47,49]]
        if x==3:
            M=[[6,50,48]]
        if x==4:
            M=[]
        if x==5:
            M=[[7,56,45]]
        if x==6:
            M=[[12,47,3]]
        if x==7:
            M=[[6,4,52],[8,37,49],[12,17,48],[18,24,25]]
        if x==8:
            M=[[9,45,57],[12,4,40]]
        if x==9:
            M=[]  
        if x==10:
            M=[[6,5,24],[12,34,16],[14,53,22],[20,3,4]]  
        if x==11:
            M=[[19,58,47]]  
        if x==12:
            M=[[8,59,48],[14,25,19]]  
           
        for i in range(0,len(M),1):   
            
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]
    
            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -20
            t2 = t1 + 80
            seis4 = Stream()
            seis4 = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis4[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc4 = seis4[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
            trc4.detrend(type='demean')
            trc4.detrend(type='linear')
            stream4.append(trc4)
            trc4.plot(type='relative',color='b')
    
    for x in range(0,len(stream4)):
        trc=stream4[x].normalize()
        stream4b.append(trc)
    
    stack_norm4 = np.sum([abs(trc.data) for trc in stream4b], axis=0)
    stack_norm4 = stack_norm4/len(stream4b)
    
    plt.figure(4)
    plt.plot(stack_norm4,color='r')
    plt.title('LS04 stacked earthquake waveform')
    
    #%% LB05
    sta = 'LS05' # STATION LB02
    cha = 'EHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(5,11,1):
        if x==5:
            M=[[7,56,46]]
        if x==6:
            M=[[12,55,2]]
        if x==7:
            M=[[6,4,53],[8,37,50],[12,17,50],[18,24,26]]
        if x==8:
            M=[[9,45,58]]
        if x==9:
            M=[[10,34,42]]  
        if x==10:
            M=[[12,34,18],[14,53,22],[20,3,3]]  
      
           
        for i in range(0,len(M),1):   
            
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]
    
            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -20
            t2 = t1 + 80
            seis5 = Stream()
            seis5 = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis5[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc5 = seis5[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
            trc5.detrend(type='demean')
            trc5.detrend(type='linear')
            stream5.append(trc5)
            trc5.plot(type='relative',color='b')
    
    
    for x in range(0,len(stream5)):
        trc=stream5[x].normalize()
        stream5b.append(trc)
    
    stack_norm5 = np.sum([abs(trc.data) for trc in stream5b], axis=0)
    stack_norm5 = stack_norm5/len(stream5b)
    plt.figure(5)
    plt.plot(stack_norm5,color='r')
    plt.title('LS05 stacked earthquake waveform')
    
    #%% LB06
    sta = 'LS06' # STATION LB02
    cha = 'EHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    year2 = 2016
    month2=6
    day2=25
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year2, month2, day2, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    
    M=[[2,28,9,0],[10,23,52,2],[21,43,19,2],[1,18,52,4],[2,5,31,4],[5,42,55,6],
       [4,10,57,7],[3,46,16,8],[18,45,31,8],[1,48,55,10],[7,36,17,12],[8,26,4,12],
       [18,30,3,15],[9,12,43,16]]
           
    for i in range(0,len(M),1):   
        
        # input time and length of waveform
        h=M[i][0]
        m=M[i][1]
        s=M[i][2]
        d=M[i][3]
        #
        t1 = t0 + d*24*60*60 + h*60*60 + m*60 + s  -20
        t2 = t1 + 80
        seis6 = Stream()
        seis6 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        seis6[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
        trc6 = seis6[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
        trc6.detrend(type='demean')
        trc6.detrend(type='linear')
        stream6.append(trc6)
        trc6.plot(type='relative',color='b')
        
    
    for x in range(0,len(stream6)):
        trc=stream6[x].normalize()
        stream6b.append(trc)
    
    stack_norm6 = np.sum([abs(trc.data) for trc in stream6b], axis=0)
    stack_norm6 = stack_norm6/len(stream6b)
    plt.figure(6)
    plt.plot(stack_norm6,color='r')
    plt.title('LS06 stacked earthquake waveform')
    
##%% LB07
#sta = 'LS07' # STATION LB02
#cha = 'HHZ' # CHANNEL - Vertical
#net = 'Z4'  # Santiaguito volcano
#loc = ''    # location, it depends mostly of which network you are in. 
#
#client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#t0 = UTCDateTime(year2, month2, day2, hour1, minute1, second1) #the format is year:day_of_the_year:month
#
#
#M=[[21,43,20,2],[1,18,53,4],[5,42,57,6],[4,10,58,7],
#   [18,45,31,8],[2,58,50,30]]#,[2,5,32,4],[18,30,4,15],[2,32,37,14]]
#       
#for i in range(0,len(M),1):   
#    
#    # input time and length of waveform
#    h=M[i][0]
#    m=M[i][1]
#    s=M[i][2]
#    d=M[i][3]
#    #
#    t1 = t0 + d*24*60*60 + h*60*60 + m*60 + s  -50
#    t2 = t1 +90 
#    seis7 = Stream()
#    seis7 = client.get_waveforms(net, sta, '', cha, t1 , t2)
#    
#    seis7[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
#    trc7 = seis7[0].slice(starttime = t1 + 30  , endtime= t2 - 10)
#    trc7.detrend(type='demean')
#    trc7.detrend(type='linear')
#    stream7.append(trc7)
##    trc7.plot(type='relative',color='b')
#
#for x in range(0,len(stream7)):
#    trc=stream7[x].normalize()
#    stream7b.append(trc)
#    
#
#stack_norm7 = np.sum([abs(trc.data) for trc in stream7b], axis=0)
#stack_norm7 = stack_norm7/len(stream7b)
#plt.figure(7)
#plt.plot(stack_norm7,color='r')
#plt.title('LB07 stacked earthquake waveform')
#%%
    return(stack_norm1,stack_norm2,stack_norm3,stack_norm4,stack_norm5,stack_norm6)#,stack_norm7)
    
    
