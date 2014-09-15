# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 13:13:03 2014

@author: ngalin
"""

#Testing of 'interpResamp_FUNC.py'
import math
import matplotlib.pyplot as plt
from interpResamp_FUNC import interpResamp

inputData = []

path = '/Users/ngalin/Documents/MATLAB/'
filename = path+'x.txt'
fin = open(filename, "r")
for line in fin:
    inputData.append(float(line.strip()))
fin.close()

tInput = list(range(0,len(inputData),1))

#test that gives back x, when I = 1, and D = 0
#I = 1; D = 0; smoothing = False
[tOutput,outputData] = interpResamp(inputData,False,1,0)
plt.plot(tInput,inputData)
plt.plot(tOutput,outputData)
if (rmse(inputData,outputData) < 1e-6) :
    print "Pass test 1"
 
#test that when interpolate to a large value, corresponding samples still overlap
#I = 120; D = 43.2; smoothing = False
[tOutput, outputData] = interpResamp(inputData,False,100,43.2)
plt.plot(tInput,inputData,'ro')
plt.plot(tOutput,outputData)
numIdx = 0
listA = []
listB = []
for t in tInput:
    try :
        idx = tOutput.index(tInput[t])
    except :
        continue
    listA.append(inputData[tInput[t]])
    listB.append(outputData[idx])
    numIdx = numIdx + 1

if (rmse(listA,listB) < 1e-6) :
    print "Pass test 2"
    
#test that when interpolate to a large value, corresponding samples still overlap
#I = 120; D = -43.2; smoothing = False
[tOutput, outputData] = interpResamp(inputData,False,100,-43.2)
plt.plot(tInput,inputData,'ro')
plt.plot(tOutput,outputData)
numIdx = 0
listA = []
listB = []
for t in tInput:
    try :
        idx = tOutput.index(tInput[t])
    except :
        continue
    listA.append(inputData[tInput[t]])
    listB.append(outputData[idx])
    numIdx = numIdx + 1

if (rmse(listA,listB) < 1e-6) :
    print "Pass test 3"


def rmse(listA,listB) :
    N = len(listA)
    M = len(listB)
    
    if (N != M):
        print "Error: RMSE received lists of different lengths"
        return 'NaN'
    else :
        diff = [(a-b)**2. for a,b in zip(listA,listB)]
        return math.sqrt(sum(diff)/N)

        