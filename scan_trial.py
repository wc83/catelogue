step=2
shift=100
crit=0.6
crite=0.6
from obspy.core import read
from obspy.signal.cross_correlation import xcorr
import glob
import time   
import obspy
   #%% constants    
shift=50*step
#%% Read in sample waveforms with which to compare all waveforms      
sample = read('/media/will/8209dcbc-09a8-4e6e-a048-6e0d2f2d2407/330a.mseed')
trace=sample[2]

#Earthquake1
start1= sample[2].stats.starttime + 3*60*60 + 52*60 + 30 #time window start - 2s before Pwave arrival
trs1 = trace.slice(starttime = start1  , endtime= start1 + 40) #cut out sample waveform with same window length as chosen event
trs1.detrend(type='demean')
trs1.filter("lowpass" ,freq=2.9)
trs1_e = obspy.signal.filter.envelope(trs1.data)

#Explosion
start2= sample[2].stats.starttime + 6*60*60 + 15*60 + 47 #time window start - 2s before Pwave arrival
trs2 = trace.slice(starttime = start2  , endtime= start2 + 40) #cut out sample waveform with same window length as chosen event
trs2.detrend(type='demean')
trs2.filter("lowpass" ,freq=2.9)
trs2_e = obspy.signal.filter.envelope(trs2.data)

#Earthquake2
start3= sample[2].stats.starttime + 6*60*60 + 1*60 + 27 #time window start - 2s before Pwave arrival
trs3 = trace.slice(starttime = start3  , endtime= start3 + 40) #cut out sample waveform with same window length as chosen event
trs3.detrend(type='demean')
trs3.filter("lowpass" ,freq=2.9)
trs3_e = obspy.signal.filter.envelope(trs3.data)

#%%
T1= time.clock()
#create empty stream for traces to be scanned
stream = sample.copy()
stream = stream.clear()
#open files to write event times into
file_quake = open('/media/will/8209dcbc-09a8-4e6e-a048-6e0d2f2d2407/catelogue/equake_test5.txt','w')
file_quake2 = open('/media/will/8209dcbc-09a8-4e6e-a048-6e0d2f2d2407/catelogue/equake2_test5.txt','w')
file_exp = open('/media/will/8209dcbc-09a8-4e6e-a048-6e0d2f2d2407/catelogue/explosion_test5.txt','w')

#%% read in day of interest and select first vertical channel, open 
for file in glob.glob('/media/will/8209dcbc-09a8-4e6e-a048-6e0d2f2d2407/*.mseed'):
    seis = read(file)
    for t in range(2,4,1):
        if seis[t].stats.channel in ['HHZ']:
            tr=seis[t]
            tr.detrend(type='demean')
            tr.filter("lowpass", freq=2.9)
            stream.append(tr)
stream.sort(keys=['starttime']) 
print('files read in') 
print(stream)
#%% cross correlation in shifting time window and find well correlated waveforms
# 
for x in range(0,len(stream),1):
    begin = stream[x].stats.starttime
    add=0          
    for t in range(0,86360, step):
        start = begin + t + add
        if t +add > 86360:
            file_exp.write("\n")
            file_quake.write("\n")
            break
        #may need to rearrange these next two lines higher into the loops to speed up CPU time
        trc = stream[x].slice(starttime = start  , endtime= start + 40)
        trc_e = obspy.signal.filter.envelope(trc.data)
        corel1=xcorr(trs1_e, trc_e, shift, full_xcorr=False)
        corel2=xcorr(trs2_e, trc_e, shift, full_xcorr=False) 
        corel3=xcorr(trs3_e, trc_e, shift, full_xcorr=False)
       
        if abs(corel1[1]) > crit:
            eventt = start - corel1[0]/100 +2
            add +=40
            file_quake.write(str(eventt))
            file_quake.write("\n")
        if abs(corel2[1]) > crite:
            expt = start -corel2[0]/100 +2
            add += 40    
            file_exp.write(str(expt))
            file_exp.write("\n")
        if abs(corel3[1]) > crit:
            expt = start -corel3[0]/100 +2
            add += 40    
            file_quake2.write(str(expt))
            file_quake2.write("\n")            
file_exp.close()
file_quake.close()
file_quake2.close()
print("End of Scan")
T2=time.clock()
elapsed= T2-T1
print('Time taken:', elapsed)
