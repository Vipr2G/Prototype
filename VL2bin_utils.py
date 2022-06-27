#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:58:35 2021

@author: russelmiller
"""

from math import trunc
from VL2shared_utils import feql, trunc2
from VL2db_utils import log_event
from VL2config import min_clstr_volume


     
def make_bins(clstr_min,clstr_max,horizon,incr):
    #note all values should be in 0.1 incr so round vs trunc?
    n_bins = round((clstr_max - clstr_min + 2 * horizon)/incr)
    bins = [0] * n_bins
    return bins
    
def init_bins(data, clstr_min, clstr_max, horizon, incr):
    bins = make_bins(clstr_min,clstr_max,horizon,incr)
    for val in data:
        bins = update_bin_cnts(val,bins,clstr_min,clstr_max, horizon, incr)
    return bins
    

def get_bin_idx(val, clstr_min, clstr_max, horizon, incr):
    # finds the bin_cnt index to update for a particular value
    val = float(val)
    c_min = float(clstr_min)
    c_max = float(clstr_max)
    horizon = float(horizon)
    incr = float(incr)
    bin_min = trunc2(c_min - horizon, incr)
    hi_bin_max = trunc2(c_max + horizon, incr)  # "hi_bin_max" is the implied upper edge, but edges are based on lo_vals only    
    if val - bin_min >= 0 and val < hi_bin_max:
        x = (val-bin_min)/incr
        if x % 1 > 0.9999:
            n = round(x)
        else:
            n = trunc(x)
        return n
    
def update_bin_cnts(val,bins,clstr_min,clstr_max,horizon,incr):
    localbins = bins.copy()
    indx = get_bin_idx(val,clstr_min,clstr_max,horizon,incr)
    #print(indx)
    if indx != None and indx<len(localbins):
        #localbins[indx] +=1
        localbins[indx] += round(1/incr)
        #print(localbins)
    return localbins    
        
def print_cnts(bins,clstr_min,clstr_max,horizon,incr):
    edges = []
    lo_val = trunc2(clstr_min - horizon,incr)
    hi_val = trunc2(lo_val + incr,incr)
    print(" bin_min to < bin_max      Count")
    for i in bins:
        edges.append([lo_val,hi_val])
        print(str(lo_val) + " to < " + str(hi_val) + "   " + str(i))
        lo_val = trunc2(lo_val + incr,incr)
        hi_val = trunc2(lo_val + incr,incr)
        

def bin_mean(bin_cnts,bin_mid_points):
    N = sum(bin_cnts)
    products = []
    for num1, num2 in zip(bin_cnts,bin_mid_points):
        products.append(num1 * num2)
    num = sum(products)
    bin_mean = num/N
    return(bin_mean)

def bin_median(bin_cnts,bin_mid_points):
    N = sum(bin_cnts)
    if N == 0:
        err1 = len(bin_cnts)
        event_str = "WARNING: Zero bin_counts VLbin_utils.bin_median len(bin_cnts) = " + str(err1)
        print(event_str)
        log_event(event_str)
        return 0,0,0
    delta = bin_mid_points[1] - bin_mid_points[0]
    median_pos = N/2
    running_count = 0
    idx = 0
    for cnt in bin_cnts:
        running_count = running_count + cnt
        if running_count >= median_pos:
            lo_cnt = median_pos - (running_count - cnt)  # median_pos - previous count
            hi_cnt = running_count-median_pos
            break
        idx +=1
    median_ratio = lo_cnt/(lo_cnt + hi_cnt)
    #median_ratio = (median_pos - lo_cnt)/(hi_cnt-lo_cnt)
    b_median = bin_mid_points[idx] + delta*(median_ratio - 0.5)
    return b_median,idx,median_ratio        
        
def bin_resize(bins,incr,old_min,new_min,old_max,new_max):
    # add or remove bins on left or right based on changes in bin_min and/or bin_max
    updated_bins = bins.copy()
    if new_min < old_min:   #add bins to left
        num_bins = round((old_min-new_min)/incr)
        add_bins = [0] * num_bins
        updated_bins = add_bins + updated_bins
    elif new_min > old_min:  #delete bins from the left
        num_bins = round((new_min-old_min)/incr)
        updated_bins = updated_bins[(num_bins):]
    if new_max > old_max:  #add bins to the right
        num_bins = round((new_max-old_max)/incr)
        add_bins = [0] * num_bins
        updated_bins = updated_bins + add_bins
    elif new_max < old_max:  # remove bins from the right
        num_bins = round((old_max-new_max)/incr)
        updated_bins = updated_bins[:(-num_bins)]
    return updated_bins
    
# calculate clstr limits from bins

def calc_bin_parms(bins,bin_min,incr,p=0.9):
    # Simplified version of calc_bin_parms...the original concept was to use the grouped mean as a measure
    # of center for the bins,but has since evolved into using the median.  In the latter case the characterization 
    # is symmetric and independent analyses of the upper/lower halves is not necessary (ASP-->SP)
    # Note: the interpretation of p is still the proportion of data on either side of center vs the entire set
    
    p_upper = .5 * (1+p)
    p_lower = 1 - p_upper
    bin_sum = sum(bins)
    
    lo_cnt = bin_sum * p_lower
    hi_cnt = bin_sum * p_upper
    
    running_cnt = 0
    lo_flag = True
    for i in range(0,len(bins)):
        running_cnt += bins[i]
        if lo_flag and lo_cnt <= running_cnt:
            lo_idx = i
            lo_flag = False
        if hi_cnt <= running_cnt:
            hi_idx = i
            break
    
    precision = 1/round(10/incr)
    clstr_min = trunc2(bin_min + lo_idx * incr,precision)
    clstr_max = trunc2(bin_min + (hi_idx + 1) * incr,precision)
    return clstr_min,clstr_max


def calc_bin_parms_old(bins,bin_min,incr,p=0.9):
    # no longer used, but preserved for posterity...if grouped mean used as the measure of center, this fctn could be
    # tweaked to achieve the original ASP intent, but as long as median is used no need to independently examine each half
    
    if len(bins) < 5:
        event_str = "INFO: unusually small bins passed to calc_bin_parms, len = "+ str(len(bins))
        log_event(event_str)
    precision = 1/round(10/incr)
    bin_mid_points = []
    delta = trunc2(incr/2,precision)
    for i in range(1,len(bins)+1):
        next_val = bin_min + (2*i-1) * delta
        bin_mid_points.append( trunc2(next_val, precision) )
    #b_mean = bin_mean(bins,bin_mid_points)
    
    if sum(bins) == 0:
        print("pause for debug (VLbinutils.calc_bin_parms)")
    
    b_median,median_idx,median_ratio = bin_median(bins,bin_mid_points)
    N_upper = sum(bins[(median_idx+1):]) + (1-median_ratio)*bins[median_idx]
    N_lower = sum(bins[0:median_idx]) + median_ratio*bins[median_idx]     #upper position not included in sum, hence median_idx vs median_idx-1
    upper_count = (1-median_ratio)*bins[median_idx]
    if upper_count >= p*N_upper:
        clstr_max_idx = median_idx
    else:
        for i in range(median_idx + 1,len(bins)):
            upper_count = upper_count + bins[i]
            if upper_count >= p * N_upper:
                clstr_max_idx = i
                break
            elif i == len(bins):   # added 210725
                clstr_max_idx = len(bins)-1   # 210727 was an error (clstr_bins vs bins) prior change didn't take??
    #print("Nupper,upper_count",N_upper,upper_count)
    c_max = trunc2(bin_mid_points[clstr_max_idx] + incr/2, incr)
    
    lower_count = median_ratio*bins[median_idx]
    clstr_min_idx = 0
    if lower_count >= p * N_lower:
        clstr_min_idx = median_idx
    else:
        for i in range(median_idx-1, 0,-1):
            lower_count = lower_count + bins[i]
            if lower_count >= p * N_lower:
                clstr_min_idx = i
                break
            elif i == 1:  # added 210725
                clstr_min_idx = 0

    c_min = trunc2(bin_mid_points[clstr_min_idx] - incr/2,incr)
    
    
    return c_min, c_max


def append_clstr_bins(data,clstrs,horizon,incr,p):
    """ 
    Given data and list [clstr_min & clstr max], create bins, reconcile clstr_min,clstr_max 
    and return appended list of [clstr_min, clstr_max, [bin_counts] ] 
    """
    clstr_inserts = []
    clstr_id = None
    
    for clstr in clstrs:
        if len(clstr) == 2:
            clstr_min = clstr[0]
            clstr_max = clstr[1]
        elif len(clstr) == 3:  #clstrs have a clstr_id as first element
            clstr_id = clstr[0]
            clstr_min = clstr[1]
            clstr_max = clstr[2]
        clstr_bin = init_bins(data, clstr_min, clstr_max, horizon, incr)
        
        # 210730 added to solve the zero bin warning....ToDo: also need to add something similar to cluster_parms
        if sum(clstr_bin) < min_clstr_volume:
            continue
            
        # compute clstr_min,clstr_max from the bin
        clstr_bin_min = trunc2(clstr_min - horizon,incr)
        
            
        new_min,new_max = calc_bin_parms(clstr_bin,clstr_bin_min,incr,p)
        if not feql(clstr_min,new_min):
            #print("changing clstr min from/to",clstr_min,new_min)
            clstr_min = new_min
        #else:
            #print("no deviation: clstr_min")
        if not feql(clstr_max,new_max):
            #print("changing clstr max from/to",clstr_max,new_max)
            clstr_max = new_max
        #else:
            #print("no deviation: clstr_max")
        
        if clstr_id == None:
            clstr_inserts.append([clstr_min,clstr_max,clstr_bin])
        else:
            clstr_inserts.append([clstr_id,clstr_min,clstr_max,clstr_bin])
    return clstr_inserts

 
    


























