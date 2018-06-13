#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 17 09:53:23 2018

@author: william
"""
def stacked_LBz_eq(fmi,fma):
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
    sta = 'LB01' # STATION LB02
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(0,13,1):
        if x==0:
            M=[[1,57,18],[8,30,14],[18,42,56],[18,59,57]]
        if x==1:
            M=[[9,16,51],[11,24,10],[13,12,8]]
        if x==2:
            M=[[3,52,28],[5,14,22],[6,1,32],[8,13,8],[9,47,48],[19,15,43]]
        if x==3:
            M=[[6,50,47],[8,4,14]]
        if x==4:
            M=[[4,26,17]]
        if x==5:
            M=[[7,56,44]]
        if x==6:
            M=[[12,47,4]]
        if x==7:
            M=[[6,4,52],[8,37,49],[12,17,49],[18,24,26]]
        if x==8:
            M=[[9,45,57],[12,4,40]]
        if x==9:
            M=[[0,4,49],[10,34,41]]  
        if x==10:
            M=[[6,5,25],[12,34,16],[14,53,20],[20,3,4]]  
        if x==11:
            M=[[19,58,46]]  
        if x==12:
            M=[[8,59,48],[14,25,20]]  
           
        for i in range(0,len(M),1):   
            
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]

            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -20
            t2 = t1 + 90
            seis1 = Stream()
            seis1 = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis1[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc1 = seis1[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
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
#    plt.title('LB01 stacked earthquake waveform')
    
    #%% LB02
    sta = 'LB02' # STATION LB02
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(0,13,1):
        if x==0:
            M=[[1,57,18],[8,30,14],[18,42,55],[18,59,58]]
        if x==1:
            M=[[10,28,22]]
        if x==2:
            M=[[3,52,28],[5,14,21],[6,1,31],[8,13,8],[9,47,48],[19,15,42]]
        if x==3:
            M=[[6,50,47],[12,23,45]]
        if x==4:
            M=[[4,26,15]]
        if x==5:
            M=[[7,56,44]]
        if x==6:
            M=[[12,47,2],[12,55,0]]
        if x==7:
            M=[[8,37,46],[12,17,47]]
        if x==8:
            M=[[9,46,3],[12,4,40]]
        if x==9:
            M=[[0,4,49],[10,34,41]]  
        if x==10:
            M=[[14,53,22],[20,3,2]]  
        if x==11:
            M=[[19,58,46]]  
        if x==12:
            M=[[8,59,48],[14,25,20]]  
           
        for i in range(0,len(M),1):   
            
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]

            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -20
            t2 = t1 + 90
            seis2 = Stream()
            seis2 = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis2[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc2 = seis2[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
            trc2.detrend(type='demean')
            trc2.detrend(type='linear')
            stream2.append(trc2)
    #        trc.plot(type='relative',color='b')
    
    for x in range(0,len(stream2)):
        trc=stream2[x].normalize()
        stream2b.append(trc)
    
    stack_norm2 = np.sum([abs(trc.data) for trc in stream2b], axis=0)
    stack_norm2 = stack_norm2/len(stream2b)
    
#    plt.figure(2)
#    plt.plot(stack_norm2,color='r')
#    plt.title('LB02 stacked earthquake waveform')
    
    #%% LB03
    sta = 'LB03' # STATION LB02
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(1,13,1):
    
        if x==1:
            M=[[5,6,46],[9,16,52],[11,24,10],[13,12,9]]
        if x==2:
            M=[[3,52,29],[5,14,23],[6,1,32],[8,13,8],[9,47,50],[20,48,28]]
        if x==3:
            M=[[16,42,43]]
        if x==4:
            M=[[8,56,14]]
        if x==5:
            M=[[7,56,44],[11,29,27],[18,29,3]]
        if x==6:
            M=[[12,47,4],[12,55,2],[15,37,52]]
        if x==7:
            M=[[8,37,49],[12,17,49],[18,24,26],[22,52,27]]
        if x==8:
            M=[[9,45,58]]
        if x==9:
            M=[[0,4,49],[10,34,41]]  
        if x==10:
            M=[[12,34,17],[14,53,22],[19,29,32],[20,3,4]]  
        if x==11:
            M=[[14,30,51],[19,58,46]]  
        if x==12:
            M=[[14,25,20]]  
           
        for i in range(0,len(M),1):   
            
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]

            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -20 
            t2 = t1 + 90  
            seis3 = Stream()
            seis3 = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis3[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc3 = seis3[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
            trc3.detrend(type='demean')
            trc3.detrend(type='linear')
            stream3.append(trc3)
    #        trc.plot(type='relative',color='b')
    
    
    for x in range(0,len(stream3)):
        trc=stream3[x].normalize()
        stream3b.append(trc)
    
    stack_norm3 = np.sum([abs(trc.data) for trc in stream3b], axis=0)
    stack_norm3 = stack_norm3/len(stream3b)
    
#    plt.figure(3)
#    plt.plot(stack_norm3,color='r')
#    plt.title('LB03 stacked earthquake waveform')
    
    #%% LB04
    sta = 'LB04' # STATION LB02
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(1,13,1):
    
        if x==1:
            M=[[9,16,51],[11,24,10],[13,12,8]]
        if x==2:
            M=[[3,52,28],[5,14,22],[6,1,32],[8,13,8],[9,47,48],[19,15,43]]
        if x==3:
            M=[[6,50,47],[8,4,14]]
        if x==4:
            M=[]
        if x==5:
            M=[[7,56,44]]
        if x==6:
            M=[[12,47,4]]
        if x==7:
            M=[[6,4,52],[8,37,49],[12,17,49],[18,24,26]]
        if x==8:
            M=[[9,45,57],[12,4,40]]
        if x==9:
            M=[[0,4,49]]  
        if x==10:
            M=[[6,5,25],[12,34,16],[14,53,22],[20,3,4]]  
        if x==11:
            M=[[19,58,46]]  
        if x==12:
            M=[[8,59,48],[14,25,20]]  
           
        for i in range(0,len(M),1):   
            
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]

            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -20
            t2 = t1 + 90
            seis4 = Stream()
            seis4 = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis4[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc4 = seis4[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
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
#    plt.title('LB04 stacked earthquake waveform')
    
    #%% LB05
    sta = 'LB05' # STATION LB02
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    for x in range(2,12,1):
        if x==2:
            M=[[3,52,27],[5,14,21],[8,13,8],[9,47,48],[19,15,43],[20,48,29]]
        if x==3:
            M=[[6,50,47]]
        if x==4:
            M=[[8,56,15]]
        if x==5:
            M=[[7,56,44]]
        if x==6:
            M=[[12,47,3],[12,55,0],[15,37,52]]
        if x==7:
            M=[[6,4,52],[8,37,50],[12,17,49],[18,24,26]]
        if x==8:
            M=[[9,45,58]]
        if x==9:
            M=[[0,4,50],[10,34,41]]  
        if x==10:
            M=[[6,5,25],[12,34,17],[14,53,22],[20,3,4]]  
        if x==11:
            M=[[14,30,53]]  
           
        for i in range(0,len(M),1):   
            
            # input time and length of waveform
            h=M[i][0]
            m=M[i][1]
            s=M[i][2]

            #
            t1 = t0 + x*24*60*60 + h*60*60 + m*60 + s  -20
            t2 = t1 + 90
            seis5 = Stream()
            seis5 = client.get_waveforms(net, sta, '', cha, t1 , t2)
            
            seis5[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
            trc5 = seis5[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
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
#    plt.title('LB05 stacked earthquake waveform')
    
    #%% LB06
    sta = 'LB06' # STATION LB02
    cha = 'HHZ' # CHANNEL - Vertical
    net = 'Z4'  # Santiaguito volcano
    loc = ''    # location, it depends mostly of which network you are in. 
    year2 = 2016
    month2=6
    day2=25
    
    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
    t0 = UTCDateTime(year2, month2, day2, hour1, minute1, second1) #the format is year:day_of_the_year:month
    
    
    M=[[2,28,10,0],[10,23,54,2],[21,43,19,2],[1,18,53,4],[2,5,33,4],[7,17,47,5],[5,42,56,6],
       [4,10,58,7],[3,46,17,8],[3,58,1,8],[18,45,34,8],[1,48,56,10],[7,36,19,12],[8,26,4,12],
       [2,32,32,14],[18,30,4,15],[9,12,43,16],[7,19,45,27],[2,58,47,30]]
           
    for i in range(0,len(M),1):   
        
        # input time and length of waveform
        h=M[i][0]
        m=M[i][1]
        s=M[i][2]
        d=M[i][3]
        #
        t1 = t0 + d*24*60*60 + h*60*60 + m*60 + s  -20
        t2 = t1 + 90
        seis6 = Stream()
        seis6 = client.get_waveforms(net, sta, '', cha, t1 , t2)
        
        seis6[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
        trc6 = seis6[0].slice(starttime = t1 + 10  , endtime= t2 - 10)
        trc6.detrend(type='demean')
        trc6.detrend(type='linear')
        stream6.append(trc6)
    #        trc6.plot(type='relative',color='b')
        
    
    for x in range(0,len(stream6)):
        trc=stream6[x].normalize()
        stream6b.append(trc)
    
    stack_norm6 = np.sum([abs(trc.data) for trc in stream6b], axis=0)
    stack_norm6 = stack_norm6/len(stream6b)
#    plt.figure(6)
#    plt.plot(stack_norm6,color='r')
#    plt.title('LB06 stacked earthquake waveform')
    
    #%% LB07
#    sta = 'LB07' # STATION LB02
#    cha = 'HHZ' # CHANNEL - Vertical
#    net = 'Z4'  # Santiaguito volcano
#    loc = ''    # location, it depends mostly of which network you are in. 
#    
#    client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
#    t0 = UTCDateTime(year2, month2, day2, hour1, minute1, second1) #the format is year:day_of_the_year:month
#    
#    
#    M=[[21,43,20,2],[1,18,53,4],[5,42,57,6],[4,10,58,7],
#       [18,45,31,8],[2,58,50,30]]#,[2,5,32,4],[18,30,4,15],[2,32,37,14]]
#           
#    for i in range(0,len(M),1):   
#        
#        # input time and length of waveform
#        h=M[i][0]
#        m=M[i][1]
#        s=M[i][2]
#        d=M[i][3]
#        #
#        t1 = t0 + d*24*60*60 + h*60*60 + m*60 + s  -50
#        t2 = t1 +90 
#        seis7 = Stream()
#        seis7 = client.get_waveforms(net, sta, '', cha, t1 , t2)
#        
#        seis7[0].filter("bandpass", freqmin=fmin,freqmax=fmax)
#        trc7 = seis7[0].slice(starttime = t1 + 30  , endtime= t2 - 10)
#        trc7.detrend(type='demean')
#        trc7.detrend(type='linear')
#        stream7.append(trc7)
#    #    trc7.plot(type='relative',color='b')
#    
#    for x in range(0,len(stream7)):
#        trc=stream7[x].normalize()
#        stream7b.append(trc)
#        
#
#    stack_norm7 = np.sum([abs(trc.data) for trc in stream7b], axis=0)
#    stack_norm7 = stack_norm7/len(stream7b)
##    plt.figure(7)
##    plt.plot(stack_norm7,color='r')
##    plt.title('LB07 stacked earthquake waveform')
    #%%
    return(stack_norm1,stack_norm2,stack_norm3,stack_norm4,stack_norm5,stack_norm6)#,stack_norm7)
    
    
