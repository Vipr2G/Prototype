#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 13:04:22 2021

@author: russelmiller
"""
import operator

def feql(x1,x2):
    if abs(x1-x2) <0.00001:
        return True
    else:
        return False
    
def rip(x,n):
    y = []
    for i in x:
        y.append(i[n])
    return y

def alt_trunc2(x,places):
    # could precompute places as: round(log(round(1/incr),10))
    # below not fully verified
    N = 10**places
    y = round(N * x)/N
    return y
    
    
def trunc2(x,incr):
    # intended use: incr = 0.1, 0.01, 0.001, etc.
    #   ...will work correctly for some other values, but not guaranteed
    #x = round(x,3)
    if incr>1:
        return round(x)
    N = round(1/incr)
    # 210729 was returning 104.7 for 104.651,.1 updates follow:
    #y = round(N * x)  # this is what it was
    
    y = N * x
    if y % 1 > 0.99999:
        y = round(y)
    else:
        y = int(y)
    #val = round(N*x)/N
    #val = trunc(y)/N
    val = y/N
    return val    

def add_unique(val,dat):  # dat is of form [[val1,count1],[val2,count2],...]
    i = 0
    for row in dat:
        if val == row[0]:
            dat[i][1] +=1
            return dat
        i +=1
    dat.append([val,1])
    return dat

def repl_nulls(dat,col,new_val):
    i = 0
    clean_dat = dat.copy()
    for row in dat:
        if row[col] == None:
            clean_dat[i,col] = new_val
        i +=1
    return clean_dat

def remove_nulls(dat):
    clean_dat = []
    for val in dat:
        if val != None:
            clean_dat.append(val)
    return clean_dat
    
            
def simple_merge(clstrs,horizon):
    """
    inputs
      - clstrs is list of [clstr_min, clstr_max], [clstr_id, clstr_min, clstr_max] or [mode_id, clstr_id, clstr_min, clstr_max]
      - horizon is in parm units
    returns 
      - merged clusters as list of [clstr_min,clstr_max
    """
    if len(clstrs) == 0:
        return []

    if len(clstrs[0]) == 2:
        cmin_idx = 0
    elif len(clstrs[0]) == 3:
        cmin_idx = 1
    elif len(clstrs[0]) == 4:
        cmin_idx = 2
    
    sorted_clstrs = sorted(clstrs, key=operator.itemgetter(cmin_idx))  #210730 fix..prev sorting on 0
    i = 0
    out_clstrs = []

    for clstr in sorted_clstrs:
        if i == 0:
            prev_min = float(clstr[cmin_idx])
            prev_max = float(clstr[cmin_idx+1])
            i +=1
        else:
            cmin = float(clstr[cmin_idx])
            cmax = float(clstr[cmin_idx+1])
            if feql(cmin, prev_max + horizon) or cmin < prev_max + horizon:
                prev_max = max(prev_max,cmax)
            else:
                out_clstrs.append([prev_min,prev_max])
                prev_min = cmin
                prev_max = cmax
    # capture last clstr
    out_clstrs.append([prev_min,prev_max])
    return out_clstrs
