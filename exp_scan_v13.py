
from obspy.core import read
from obspy.signal.cross_correlation import correlate
import glob
import time
import obspy
import numpy as np
from numpy import argmax
import matplotlib.pyplot as plt
   #%% constants    
shift=100
step=2
crite=0.55 # REDUCED FROM V10
window=45
fmin=0.1
fmax=10

#%% Read in sample waveforms with which to compare all waveforms      
sample = read('/Users/william/Documents/lb01/14_330z.mseed')
trace=sample[0]
trace.detrend(type='demean')
trace.filter("bandpass", freqmin=fmin,freqmax=fmax)

trs_e=stacked_exp(fmin,fmaz)
#%%
T1= time.clock()
#create empty stream for traces to be scanned
stream = sample.copy()
stream = stream.clear()
#open files to write event times into
file_exp = open('/Users/william/Documents/scanner/exp_scan_explosion_v13.txt','w')
#%% read in day of interest and select first vertical channel, open 
#nu=0
for file in glob.glob('/Users/william/Documents/lb01/14_33*.mseed'):
#    nu+=1
#    if nu==3:
#        break
    seis = read(file)
    tr =seis[0]
    tr.detrend(type='demean')
    tr.filter("bandpass", freqmin=fmin,freqmax=fmax)
    stream.append(tr)
stream.sort(keys=['starttime']) 
print('files read in') 
#print(stream)

#%% cross correlation in shifting time window and find well correlated waveforms
# 
Times_EXP=[]
week=[]
Day=[]
epd=[]
epw=[]
w=0
wn=0
event_number=0
event_num2=0
Event_store=[]
for x in range(0,len(stream),1): 
    begin = stream[x].stats.starttime
    add=0
    lastcorr=0         
    dl=stream[x].stats.endtime - stream[x].stats.starttime
        
    for t in range(0,86360, step):
        view_start = begin + t + add
        view_end = view_start + window
        if view_start > stream[x].stats.endtime - window :
#            print('new day')
            file_exp.write("\n")
#            file_exp_near.write("\n")
#            file_unk.write("\n")
            print('Day',x+1,'of',len(stream))
            break

        trc = stream[x].slice(starttime = view_start  , endtime= view_end)
        trc.detrend(type='demean')
        # only look for things above noise level 
        amp=abs(trc.max())
        
        if amp > 1600:
#            print('amp=',amp)
            #frequency info
            peak,cf,bwid,bwid25 = freq_info(trc.data,view_start,view_end)
#            tr_data=trc.data

                   
            if 0.75 < cf < 2.75:    # Changed FROM V10
#                print('cf=',cf)
                if 0.2 < peak < 2.5: # Changed FROM V10
#                    print('peak=',peak)
                    if 0.2 < bwid < 3: # Changed FROM V10
#                        print(view_start)
                        trc_e = obspy.signal.filter.envelope(trc.data) 
            #correlate between data and earthquake envelopes
                        top_v,top,corell = corel(trs_e,trc_e,shift)
                        
                        if 0.0 in corell : # in the case of missing data - which will cause function to crash
                            break            
                        
                        
            #            if the correlation is positive, exp found
                        if abs(top_v) > crite:
                            expt = view_start - (top/100) +2
                            file_exp.write(str(expt))
                            Times_EXP.append(expt)
                            file_exp.write("\n")
                            add += 90-step               #skip 90s to avoid double catch
                            lastcorr=0
                            event_number+=1
                            event_num2+=1
                            
                            rt=expt.timestamp	
                            jd=expt.julday
                            yr=expt.year
                            mo=expt.month
                            da=expt.day
                            hr=expt.hour
                            mi=expt.minute
                            se=expt.second
                            ms=expt.microsecond
                            row=([rt,jd,yr,mo,da,hr,mi,se,ms,dl])
                            Event_store.append(row)
                            
#                            print('match at ', expt)
#                            stream[x].plot(type='relative',color='b', starttime=view_start , endtime=view_end)
 
      
        else:
            add += 38
            lastcorr=0       
    epd.append(event_number)
    event_number=0
    t=x+1
    Day.append(t)
    w+=1
    if w==7:
        wn+=1
        week.append(wn)
        epw.append(event_num2)
        w=0
        event_num2=0
        
        
#%%   save files with info
print("End of Scan")
T2=time.clock()
elapsed= T2-T1
print('Time taken:', elapsed)
file_exp.close()
np.savetxt("EXP_times_v13.csv", Event_store,delimiter=",",header="Time_stamp,Day_numver,Year,Month,Day,Hour,Min,Sec,Milisec,record_length")
#np.savetxt("EXP_epd_v10.csv", epd,delimiter=",",header="day_number,events,day_length")

#%% Plot info
diff=[]
for x in range(1,len(Times_EXP)):
    t_diff=Times_EXP[x]-Times_EXP[x-1]
    diff.append(t_diff)
    diff.sort()
    diff = [x for x in diff if x <= (6*60*60)] #longer than 6hrs, considered to have missed an event    
plt.hist(diff,bins=30)
#plt.savefig('time_gap_exp_events_v11')


fig2 = plt.figure()
plt.plot(Day,epd)
plt.ylim([0,max(epd)+10])
#plt.savefig('Daily_explosions_v11')

fig3 = plt.figure()
plt.plot(week,epw)
plt.ylim([0,max(epw)+100])
#plt.savefig('Weekly_explosions_v11')


