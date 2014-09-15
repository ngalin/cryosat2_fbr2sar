"""
Created on Mon Sep 15 13:05:39 2014

@author: ngalin
"""

# -*- coding: utf-8 -*-


#import math
import cmath
import numpy as np

def interpResamp(inputData,smoothing,I,D) :

    outputData = []
    N = len(inputData)
    n = list(range(0,N,1))
    
    t = list(frange(0,N,1./I))

    t = [x+D for x in t]
    L = len(t)
    
    for i in range(0, L):
        idx = [(t[i]-x) for x in n]
        snc = (np.sinc(idx)).tolist()
        tmp = [a*b for a,b in zip(inputData,snc)]
        
        if (smoothing) :
            win = getWindowExpression(idx)
            tmp = [a*b for a,b in zip(tmp,win)]
        
        outputData.append(sum(tmp))
        
    outputData = [abs(x) for x in outputData]
    return [t, outputData]
 
#from: http://stackoverflow.com/questions/477486/python-decimal-range-step-value/20549652#20549652
def frange(start,end,step):
    return map(lambda x: x*step, range(int(start*1./step),int(end*1./step)))
    
def getWindowExpression(idx) :
    #factor of 0.4 is an empirical value - change it as you see fit
    tmp = [(complex(x)**(0.4)) for x in idx]
    return [cmath.exp(-x) for x in tmp] 
    
