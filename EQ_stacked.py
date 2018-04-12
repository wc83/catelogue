
import numpy as np
from obspy.core import read
from obspy.signal.cross_correlation import correlate
import obspy
import matplotlib.pyplot as plt

#%%

shift=300


seisb = read('/Users/william/Documents/lb01/14_326z.mseed')
stream=seisb.copy()
stream=seisb.clear()
stream2=stream.copy()
#%%
for x in range(0,14,1):
    if x==0:
        seis = read('/Users/william/Documents/lb01/14_326z.mseed')
        M=[[13,48,5],[15,3,20],[15,13,21]]
    if x==1:
        seis = read('/Users/william/Documents/lb01/14_328z.mseed')
        M=[[1,57,18],[8,30,14],[18,42,56],[18,59,58]]
    if x==2:
        seis = read('/Users/william/Documents/lb01/14_329z.mseed')
        M=[[5,6,46],[9,16,52],[10,28,22],[11,24,10],[13,12,9]]
    if x==3:
        seis = read('/Users/william/Documents/lb01/14_330z.mseed')
        M=[[3,52,27],[5,14,21],[6,1,32],[8,13,8],[9,47,48],[19,15,43],[20,48,29]]
    if x==4:
        seis = read('/Users/william/Documents/lb01/14_331z.mseed')
        M=[[6,50,47],[8,4,14],[12,23,35],[16,42,43]]
    if x==5:
        seis = read('/Users/william/Documents/lb01/14_332z.mseed')
        M=[[4,26,17],[8,56,15]]
    if x==6:
        seis = read('/Users/william/Documents/lb01/14_333z.mseed')
        M=[[7,56,44],[11,29,26],[18,29,3]]
    if x==7:
        seis = read('/Users/william/Documents/lb01/14_334z.mseed') 
        M=[[12,47,4],[12,55,0],[15,37,52]]
    if x==8:
        seis = read('/Users/william/Documents/lb01/14_335z.mseed')
        M=[[2,15,37],[6,4,52],[8,37,49],[12,17,49],[18,24,26],[22,52,27]]
    if x==9:
        seis = read('/Users/william/Documents/lb01/14_336z.mseed')
        M=[[9,45,58],[12,4,40]]
    if x==10:
        seis = read('/Users/william/Documents/lb01/14_337z.mseed')
        M=[[0,4,49],[10,34,41]]  
    if x==11:
        seis = read('/Users/william/Documents/lb01/14_338z.mseed')
        M=[[6,5,25],[12,34,17],[14,53,20],[19,29,32],[20,3,4]]  
    if x==12:
        seis = read('/Users/william/Documents/lb01/14_339z.mseed')
        M=[[14,30,53],[19,58,46]]  
    if x==13:
        seis = read('/Users/william/Documents/lb01/14_340z.mseed')
        M=[[8,59,48],[14,25,20]]  
       
    for i in range(0,len(M),1):   
        
        # input time and length of waveform
        h=M[i][0]
        m=M[i][1]
        s=M[i][2]
        window=80
        #
        seis[0].filter("bandpass", freqmin=1,freqmax=12)
        start=seis[0].stats.starttime + h*60*60 + m*60 + s -25 -30 # define start time (hours mins seconds) 
        end = start + window +60 # define length of window in seconds
        trc = seis[0].slice(starttime = start  , endtime= end)   
        trc.detrend(type='demean')
        stream.append(trc)
        
#        trc.plot(type='relative',color='b',starttime = start  , endtime= end)
#print(stream)
#%%        

plt.figure(1)
stack = np.sum([trc.data for trc in stream], axis=0)
#print(stack)
stack = stack/50
plt.plot(stack,color='g')

for x in range(0,50):
    tr=stream[x].normalize()
    stream2.append(tr)
plt.figure(2)
stack_norm = np.sum([tr.data for tr in stream2], axis=0)
stack_norm = stack_norm/50
#print(stack_norm)
plt.plot(stack_norm,color='r')

bottom=1
for p in range(0,50):
    tr = stream[p].normalize()
    corel=correlate(stack_norm, tr, shift, demean=True,normalize='naive',domain='time' )
    top=corel.argmax()
    bot=corel.argmin()
    top_v = corel[top]
    bot_v = corel[bot]  
    
    if abs(bot_v) > abs(top_v):
        top_v=abs(bot_v)
        top=100-bot #goes from -100 to +100 but list goes from 0 to 200
        #print('data-eq:', top_v)
#    print('Xc',top_v)
    if top_v < bottom:
        bottom=top_v
print('Bot Xc',bottom)


