#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 08:12:26 2021

@author: russelmiller
"""

import numpy as np
from statistics import median
from VL2shared_utils import feql, trunc2, add_unique, rip
from VL2config import min_clstr_ints, min_clstr_peak, min_clstr_volume, min_clstr_density


def get_clstr_id(val,clstrs):  # may be a duplicate fctn
    for clstr in clstrs:
        if val >= clstr[1] and val < clstr[2]:
            return clstr[0]
    return 0

def wobulate_local_seed(data,clstr_min,clstr_max,incr,horizon,p):  # newversion processes data already in memory vs going back to the database
    """ 
    documentation tbd
    """
    clstr_min = trunc2(clstr_min,incr) # shouldn't be necessary since clstrs stored in DB with ptoper precision
    clstr_max = trunc2(clstr_max,incr) # shouldn't be necessary since clstrs stored in DB with ptoper precision
 
    stable = 0
    while stable != 1:
        dat = []
        for val in data:
            if val >= clstr_min - horizon and val < clstr_max + horizon:
                dat.append(val)
           
        #dat = (val for val in data if val >= clstr_min - horizon and val < clstr_max + horizon)
        dat.sort()  # just in case...
        #print(dat)
        #stable = 1

        if dat == None or dat == []:
            print("***** WARNING: no data to support cluster (wobulate_seed2)")
            return clstr_min,clstr_max

        dat = np.array(dat).astype(float)
                #print("elnot,parm,mod_type,clstr_min,clstr_max,horizon,dat",elnot,parm,mod_type,clstr_min,clstr_max,horizon)
        #print(dat)
        c_mean = median(dat)
        #print('c_mean',c_mean)
        up_dat = []
        for x in dat:
            if x >= c_mean:
                up_dat.append(x)
        c_max = np.percentile(up_dat,100*p)
        #print('c_max',c_max)
        dn_dat=[]
        for x in dat:
            if x <= c_mean:
                dn_dat.append(x)
        c_min = np.percentile(dn_dat,100*(1-p))
        c_min,c_max = trunc2(c_min,incr), trunc2(c_max+incr,incr)
        #print('c_min',c_min)
        #if c_min == clstr_min and c_max == clstr_max:
        if feql(c_min,clstr_min) and feql(c_max,clstr_max):
            stable = 1 # prob not needed
            return c_min,c_max
        else:
            clstr_min = c_min
            clstr_max = c_max



def sim_seq(sq1,sq2,mod_type):
    # compare 2 sequences with same mod_type
    # input format sq1,sq2 both consist of
    if len(sq1) == 0 or len(sq2) == 0:
        return False

    if mod_type == 'D':
        if sq1 == sq2:
            return True
        else:
            return False
    elif mod_type == 'M':
        # require child to be a subsequence of parent
        if len(sq1) > len(sq2):
            parent = sq1.copy()
            child = sq2.copy()
        else:
            parent = sq2.copy()
            child = sq1.copy()
        indices = []
        for i in range(len(parent)):
            if child[0] == parent[i]:
                indices.append(i)
        for index in indices:
            is_found = True
            parent_index = index
            for item in child:
                if item != parent[parent_index]:
                    is_found = False
                parent_index += 1
                if parent_index >= len(parent):
                    parent_index = 0
            if is_found:
                return True
        return False
    elif mod_type == 'S' or mod_type == 'J':
        #require all elements of child to be in parent
        if len(sq1) > len(sq2):
            parent = sq1.copy()
            child = sq2.copy()
        else:
            parent = sq2.copy()
            child = sq1.copy()
        for elem in child:
            for vals in parent:
                if elem == vals:
                    match = 1
                    break
                else:
                    match = 0
            if match == 0:
                return False
        return True


def encode_int_pris(dat, pri_clstrs):
    if len(dat) == 0:
        print("***** No PRI data for sequence building")
        return []
    if len(pri_clstrs) == 0:
        print("***** No PRI clusters for sequence building")
        return []
    # Represent intercept sequences in terms of clstr_id
    int_seqs = []
    cur_seq = []
    bad_seq = False
    for row in dat:
        int_id = row[0]
        pval = float(row[1])
        pnum = row[2]
        
        clstr_id = get_clstr_id(pval,pri_clstrs)
        #clstr_id = get_pri_clstr(elnot,mt,pval)

        if pnum == 1:   # --> starting a new sequence
                
            if cur_seq != []:  # cur_seq is a prior sequence that may need to be captured
                # save prev seq if it exists and is valid
                if bad_seq == False:
                    int_seqs.append(cur_seq)
            # back to dealing with current sequence, (pnum = 1)
            bad_seq = False
            if clstr_id == 0:   # evaluating first element (pnum = 1)
                bad_seq = True
            cur_seq = [int_id,clstr_id]   #starting the new sequence
        else:  #pnum != 1
            if clstr_id == 0:
                bad_seq = True
            cur_seq.append(clstr_id)
    #capture last sequence
    if bad_seq == False:
        int_seqs.append(cur_seq)

    return int_seqs


def get_residue_seqs(dat,pri_clstrs,min_seq_cnt):
    # PURPOSE: convert intercept PRI sequences to sequences of clstr_id (coded pri sequences)
    #          return only unique coded sequences with intercept associations
    
    int_seqs = encode_int_pris(dat,pri_clstrs)
       
    # Find unique sequences meeting threshold
    x2 = [[0],[0,0]]  #np.unique runs amok if all seqs are same length
    for i in range(0,len(int_seqs)):
        x2.append(int_seqs[i][1:])
    y = np.array(x2, dtype=object)
    #print("y",y)
    z = np.unique(y)  # returns array of lists unless all lists are of length 1, then returns array of int
    #print("z",z)
    out_seqs = []
    seq_num = 1
    for seqs in z:
        seq_cnt = sum(1 for i in y if i == seqs)  # note can alternatively return counts from np.unique()
        if seq_cnt >= min_seq_cnt and seqs != [0] and seqs != [0,0]:
            out_seqs.append([seq_num,seqs,seq_cnt])
            seq_num +=1
    # get int-seq associations
    seq_ints = []
    for int_seq in int_seqs:
        for seq in out_seqs:
            if int_seq[1:] == seq[1]:
                #seq_ints.append([int_seq[0],seq[0]])  #returns [intercept_id, seq_id]
                seq_ints = add_unique([int_seq[0],seq[0]],seq_ints)   # modified 210725 to switch to add_unique....didn't solve, orig prob ok
    seq_ints = rip(seq_ints,0)
    return out_seqs,seq_ints
    ####################
    # END get_residue_seqs
    ####################          

def group_seqs(seqs,mod_type):   # seqs = list of of [seq_num, [seq], seq_cnt]
    max_len = 0
    for seq in seqs:
        l = len(seq[1])
        if l > max_len:
            max_len = l


    sorted_seqs = []
    for length in range(max_len,0,-1):
        for seq in seqs:
            l = len(seq[1])
            if l == length:                
                sorted_seqs.append([ seq[0],seq[1],seq[2]])
    #print("sorted_seqs",sorted_seqs)
    # now group similar sequences
    parent_seqs = []
    child_seqs = []
    i = 0
    for seq in sorted_seqs:
        i +=1
        sq1 = seq[1]
        #sq1_cnt = seq[2]
        skip_seq = False
        for child in child_seqs:
            if seq[1] == child[1]:
                # seq already disposed of
                skip_seq = True
                break
        if skip_seq == False:
            parent_seqs.append(seq)
            compare_seqs = sorted_seqs[i:]
            for seq2 in compare_seqs:
                sq2 = seq2[1]
                sq2_cnt = seq2[2]
                if sim_seq(sq1,sq2,mod_type) == True: # and sq1 != sq2: # 210703 @1639 -- prevent creating identical children
                    child_seqs.append([seq[0],sq2,sq2_cnt,seq2[0]])  #[parent_id, [child_seq], child_hc, add: c_id for int association]

    return parent_seqs,child_seqs
##### End group_seqs()


def cluster_data(data,incr,horizon,thresh):
    if data == [] or data == None or len(data) < min_clstr_ints:
        return []
    data.sort()
    parm_dat = np.array(data).astype(float)
    
    # brute force parm slope (no processing of gaps, appears more efficient than scipy.stats.cumfreq)
    x_min = min(parm_dat)
    #print(x_min)
    
    start_val = trunc2(x_min,incr)
    
    n = 0
    parm_slope=[]
    max_slope = 0.0
    sum_n = 0 #testing
    for val in parm_dat:
        if val < (start_val + incr):
            n +=1
        else:
            slope = n/incr
            if slope > max_slope:
                max_slope = slope
            parm_slope.append([start_val,slope])
            sum_n = sum_n + n
            n = 1
            start_val = trunc2(val,incr)
    # capture last slope
    slope = n/incr
    if slope > max_slope:
        max_slope = slope
    parm_slope.append([start_val,slope])
    sum_n = sum_n + n
    #print('sum_n',sum_n)
    #print('parm_slope,parm,mt',parm,mod_type,parm_slope)    
    
    if max_slope < min_clstr_peak/incr:
        return []
    
    #print("max_slope: ",max_slope)
    #print(len(data),len(parm_slope))
    #print(parm_slope)
    # Extract clusters from parm_slope
    out_clstrs = []
    in_clstr = False
    in_lull = False
    avg_slope = 0
    k = 0
    last_valid = 0
    clstr_min = -1
    clstr_count = 0
    clstr_count2 = 0
    csmax = 0  #220202 "clstr slope max" all references added 220202
    
    for val in parm_slope:
        if not in_clstr:
            if val[1] < thresh * max_slope:
                continue
            else:
                in_clstr = True
                clstr_min = val[0]
                last_valid = clstr_min
                clstr_count = val[1]
                if val[1] > csmax:
                    csmax = val[1]
                clstr_count2 = 0
                in_lull = False
                avg_slope = 0
                k = 0
        else:  # in a clstr
            if val[0] - last_valid > horizon:   #cao of prior clstr and start a new one if thresh exceeded or continue search
                clstr_max = last_valid + incr

                ### 220202 update:
                cvolume = (clstr_count-clstr_count2)*incr
                cdensity = cvolume/(clstr_max-clstr_min)
                cpeak = csmax * incr
                if cvolume >= min_clstr_volume and cpeak >= min_clstr_peak and cdensity >= min_clstr_density:
                #if clstr_count-clstr_count2 >= min_clstr_volume/incr:
                    out_clstrs.append([trunc2(clstr_min,incr),trunc2(clstr_max,incr)])
                clstr_count = 0  #Bug fix 220127
                csmax = 0 # Added 220202
                
                if val[1] >= thresh * max_slope:
                    clstr_min = val[0]
                    last_valid = val[0]
                    clstr_count += val[1]
                    if val[1] > csmax:
                        csmax = val[1]
                    clstr_count2 = 0
                    avg_slope = 0
                    in_clstr = True
                    in_lull = False
                    continue
                else:  
                    clstr_count = 0
                    clstr_count2 = 0
                    in_clstr = False
                    in_lull = False
                    avg_slope = 0
                    k = 0
                    continue
            #prior continues make this an else: still within horizon
            if not in_lull:
                if val[1] >= thresh * max_slope:
                    last_valid = val[0]
                    clstr_count += val[1]
                    if val[1] > csmax:
                        csmax = val[1]
                    avg_slope = 0
                    continue
                else: #not in a lull yet but lost threshold
                    clstr_count += val[1]
                    if val[1] > csmax:
                        csmax = val[1]
                    clstr_count2 = val[1]
                    k = round((val[0]-last_valid)/incr)
                    avg_slope = val[1]/k
                    #k = 1
                    in_lull = True
                    continue
            else:  # already in a lull
                clstr_count += val[1]
                if val[1] > csmax:
                    csmax = val[1]    # in this case, could result in higher peak than actually in clstr, but should be rare 
                clstr_count2 += val[1]
                k = round((val[0]-last_valid)/incr)
                avg_slope = ( (k-1)*avg_slope + val[1] ) / k  # Updated 220127: k-1/k vs k/k+1 since k already updated
                                                              # Note: empty incrs get credit for prior avg_slope

                if avg_slope >= thresh * max_slope:  
                    clstr_count2 = 0
                    last_valid = val[0]
                    in_lull = False
                    avg_slope = 0
                    k = 0
    cvolume = (clstr_count - clstr_count2) * incr
    cdensity = cvolume/(last_valid+incr - clstr_min)
    cpeak = csmax * incr     
    if in_clstr and cvolume >= min_clstr_volume and cpeak >= min_clstr_peak and cdensity >= min_clstr_density:           
    #if in_clstr == True and clstr_count-clstr_count2 >= min_clstr_volume/incr:
        clstr_max = last_valid + incr
        out_clstrs.append([trunc2(clstr_min,incr),trunc2(clstr_max,incr)])

    return out_clstrs

def purge_conformal(data,clstrs,horizon):
    """
    Given a list of data values and a list of clusters removes data points falling within a horizon of any cluster
    
    # 211006: commented out horizon buffer to prevent fighting between residue miner and associator
    
       - clstrs is a list of [clstr_id, clstr_min, clstr_max] or just [clstr_min, clstr_max]
       - horizon is absolute and in parameter units
    Returns the list of data points that did NOT match any cluster
    """
    ff = 1e-5  #fudge factor to prevent boundary values from lingering in residue
    residue = []
    for row in data:
        if type(row) == tuple:
            value = float(row[0])
        else:
            value = float(row)
        match = False
        for clstr in clstrs:
            if len(clstr) == 3:
                c_min = float(clstr[1])  # - horizon
                c_max = float(clstr[2])  # + horizon
            elif len(clstr) == 2:
                c_min = float(clstr[0])  # - horizon
                c_max = float(clstr[1])  # + horizon
            if value+ff >= c_min and value-ff < c_max:
                match = True
                break
        if match == False:
            residue.append(value)
    return residue


def purge_conformal2(data,clstrs,horizon):
    """
    Same as purge_conformal, but input includes and process retains int_id's
    
    Given a list of data values and a list of clusters removes data points falling within a horizon of any cluster
       - clstrs is a list of [clstr_id, clstr_min, clstr_max] or just [clstr_min, clstr_max]
       - horizon is absolute and in parameter units
    Returns the list of data points that did NOT match any cluster
    """
    if len(clstrs) == 0 or len(data) == 0:
        return rip(data,0)
    
    residue = []
    res_ids = []
    if len(clstrs[0]) == 2:
        cmin_idx = 0
    else:
        cmin_idx = 1
    for row in data:
        int_id = row[0]
        value = row[1]
        match = False
        for clstr in clstrs:
            c_min = float(clstr[cmin_idx]) #- horizon
            c_max = float(clstr[cmin_idx+1]) #+ horizon
            if value >= c_min and value < c_max:
                match = True
                break
        if match == False:
            residue.append(value)
            res_ids = add_unique(int_id,res_ids)
    res_ids = rip(res_ids,0)
    return res_ids

def get_conformal(data,clstrs,horizon):
    """
    Opposite of purge_conformal, returns int_id's for intercept values that fit inside the provided clusters
    
    Given a list of data values and a list of clusters removes data points falling within a horizon of any cluster
       - clstrs is a list of [clstr_id, clstr_min, clstr_max] or just [clstr_min, clstr_max]
       - horizon is absolute and in parameter units
    Returns the list of data points that did NOT match any cluster
    """
    if len(clstrs) == 0 or len(data) == 0:
        #return data
        return []    # 210915 nothing should be conformal in this case?
    
    conformal_ids = []
    if len(clstrs[0]) == 2:
        cmin_idx = 0
    else:
        cmin_idx = 1
    for row in data:
        int_id = row[0]
        value = row[1]
        for clstr in clstrs:
            c_min = float(clstr[cmin_idx]) - horizon
            c_max = float(clstr[cmin_idx+1]) + horizon
            if value >= c_min and value < c_max:
                conformal_ids.append(int_id)
                break
    return conformal_ids


def true_heard_count(int_stats):   # instats is array of [collector, int_up_time (day/hour only), ...?]
    # need to ensure int_up_time provided or converted to reasonable level of granularity
    unique_collectors = []
    unique_times = []
    for row in int_stats:
        collector = row[0]
        int_up_time = row[1]
        unique_collectors = add_unique(collector)
        unique_times = add_unique(int_up_time)
    x = 0
    y = 0
    for i in unique_collectors:
        x += i[1]
    for i in unique_times:
        y += i[1]
    true_count = max(x,y)
    return true_count
