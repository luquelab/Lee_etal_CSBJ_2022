# -*- coding: utf-8 -*-
"""
Author:
    Diana Y. Lee, Luque Lab, SDSU
    dlee@sdsu.edu
Purpose:
    T-number functions reused across all G2T/MCP2T modules
"""

# basic imports
import numpy as np

def T_power(x,y):
    """
    This function takes the log of two arrays of data (power transforms both independent and dependent variables)
    input: two arrays, x and y
    output: two arrays, log(x) and log(y)
    """
    x_power = np.log(x)
    y_power = np.log(y)
    return x_power, y_power

def tModel():
    """
    This function returns the current tModel linear regression coefficients.
    input: none
   
    output: a 1x4 array of the slope, intercept, slope error and intercept error  
    
    """ 
    tModCurr = [0.71267, 0.03334, -0.71778, 0.14706]
    return tModCurr

from bisect import bisect_left

def tNearest(tPossible, tRaw):
    """
    This function finds the closest valid T-number from a list
    
    input:  tPossible - a sorted list of T-numbers (integers)
            tRaw - a decimal number to sort
    output: the closest T-number from tPossible to tRaw. If two numbers are equally close, the smaller T will be returned.   
    Assumes tPossible is sorted. Returns closest value to tRaw.
    If two numbers are equally close, the smaller T will be returned.
    """
    pos = bisect_left(tPossible, tRaw)
    if pos == 0:
        return tPossible[0]
    if pos == len(tPossible):
        return tPossible[-1]
    before = tPossible[pos - 1]
    after = tPossible[pos]
    if (after - tRaw) < (tRaw - before):
       return after
    else:
       return before

def tNearestFloor(tPossible, tRaw):
    """
    This function finds a valid T-number from a list by rounding down the raw T input 
    
    input:  tPossible - a sorted list of T-numbers (integers)
            tRaw - a decimal number to sort
    output: the closest T-number from tPossible to tRaw, rounded down 
    Assumes tPossible is sorted. 
    """
    pos = bisect_left(tPossible, tRaw)
    if pos == 0:
        return tPossible[0]
    if pos == len(tPossible):
        return tPossible[-1]
    if ((pos > 0) and (pos < len(tPossible))):
       return tPossible[pos - 1]
   
def tNearestValid(tPossible, tRaw, CI=0.15):
    """
    This function finds the closest valid T-number from a list within a given percentage confidence interval
    input:  tPossible - a sorted list of T-numbers (integers)
            tRaw - a decimal number to sort
            CI - a decimal value between 0 and 1 (percent), If none is provided, the default is 0.15
    output: the closest T-number from tPossible to tRaw if it is within the CI percentage. 
            Otherwise, returns 0. If two numbers are equally close, the smaller T will be returned.   
    """
    pos = bisect_left(tPossible, tRaw)
    if pos == 0:
        return tPossible[0]
    if pos == len(tPossible):
        before = tPossible[pos - 1]
        if tRaw <(before+before*(CI)):
            return before
        else:
            return 0
    before = tPossible[pos - 1]
    after = tPossible[pos]
    if (after - tRaw) < (tRaw - before):
        if tRaw >(after-after*(CI)):
            return after
        else:
            return 0
    else:
        if tRaw <(before+before*(CI)):
            return before
        else:
            return 0

def tNum(genomesize,tType=0, CI=0.15):
    """
    This function calculates the T-number based on the genome size.
    input: genomesize - a decimal value expressed in kilo base pairs
           tType - the T-number type desired. Valid types are 0 (raw), 1 (nearest), or 
                   2 (nearest within a confidence interval). Defaults to 0.
           CI  - a decimal value between 0 and 1 (percent), If none is provided, the default is 0.15
   
    output: T-number  
    
    """ 
    
    tType = tType or 0
    if (tType > 2):
        print("Valid tTypes are 0 (raw), 1 (round) or 2 (nearest T within a confidence interval). Defaults to 0.")
        return

    tps2, tps, tps_t, tps_h = tList(7)
    tMod = tModel()

    tRaw = np.exp(tMod[0]*np.log(genomesize)+tMod[2])
    if (tType == 0):
        return tRaw
    if (tType == 1):
        tRound = tNearest(tps,tRaw)
        return tRound
    if (tType == 2):
        tReal = tNearestValid(tps,tRaw, CI)
        return tReal

def tList(hkLim1):
    """
    This function creates possible T-number lists. T-number is equal to h**2 + hk + k**2
    input: hkLim1 - integer limit of h and k
   
    output: t1 - an array of all possible t numbers 
            t2 - an array of all possible t numbers except 0
            th - an array of all t_hex t numbers
            tt - an array of all t_trihex t numbers
    
    """ 
    hTest1 = np.zeros(hkLim1)
    kTest1 = np.zeros(hkLim1)

    for i in range (hkLim1):
        hTest1[i]=i;
        kTest1[i]=i+1;

    tps1 = [];
    tps_t1 = [];
    tps_h1 = [];

    for i in range(hkLim1):
        for j in range(hkLim1):
            tps1.append(i**2+i*j+j**2)
            tps1.append(round((i**2+i*j+j**2)*(4/3),2))
            tps_t1.append(i**2+i*j+j**2)
            tps_h1.append(round((i**2+i*j+j**2)*(4/3),2))

    t1 = np.unique(np.asarray(tps1));
    t2 = np.unique(np.asarray(tps1))[1:];
    tt = np.unique(np.asarray(tps_t1))[1:];
    th = np.unique(np.asarray(tps_h1))[1:];
    
    return [t1, t2, tt, th]

def tDictAll(hkLim2):
    """
    This function creates dictionaries for possible T-number lists with integer values
    input: hkLim2 - integer limit of h and k
   
    output: tdictionary - dictionary that yields integer indices for each possible T-number 
            tdictrev - dictionary that yields the T-number associated with integer index
    
    """ 
    #create a list of possible, valid T-numbers, as well as separate t-number lists for T_h and T_t through h/k = hkLim2 
    tps2, tps, tps_t, tps_h = tList(hkLim2)

    # create a numbered list of the posible T-numbers
    tIndex = np.arange(len(tps))
    tIndex

    # create a zip object from the two lists above
    zipbObj = zip(tps, tIndex)

    # create a dictionary from zip object
    tdictionary = dict(zipbObj)
    tdictionary.update( {0:len(tps)} )

    # and reverse it
    zipbObj = zip(tIndex, tps)
    tdictrev = dict(zipbObj)
    tdictrev.update( {len(tps):0} )
    
    return [tdictionary,tdictrev]
