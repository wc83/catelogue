#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 10:08:56 2018

@author: william
"""

def corel(wav1,wav2,shift):
    from obspy.signal.cross_correlation import correlate

    corel=correlate(wav1, wav2, shift, demean=True,normalize='naive',domain='time' )
    top=corel.argmax()
    bot=corel.argmin()
    top_v = corel[top]
    bot_v = corel[bot]  
    
    if abs(bot_v) > abs(top_v):
        top_v=abs(bot_v)
        top=shift-bot #goes from -shift to +shift but list goes from 0 to 2xshift

    return(top_v,top,corel)