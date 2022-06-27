#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 09:11:42 2021

@author: russelmiller
"""
from VL2db_utils import db_con, get_tuning, db_insert_seqs,db_insert_local_parms
from VL2config import min_clstr_peak, min_clstr_volume, fq_idb,min_residue_ints,min_seq_cnt
from VL2shared_utils import add_unique, rip, trunc2, simple_merge, feql
from VL2dm_utils import get_residue_seqs, group_seqs
from VL2associator import new_associator
#import numpy as np
#min_res_ints = 1
k_buf = 1  #horizon multiplier for bin buffer...integrate into VLconfig and leverage everywhere necessary


def new_residue_miner(notation = None, mine_clstrs = True):
    # 1. Mine global clusters
    # 2. Form/extract new modes
    
    db = db_con()
    cur = db.cursor()
    mode_parm_inserts = []
    
    ## 1. Mine & store new clusters
    if mine_clstrs:   # can skip cluster mining when called by baseline_builder
        if notation == None:
            cur.execute("select distinct elnot, mod_type from parm_residue")
        else:
            cur.execute("select distinct elnot,mod_type from parm_residue where elnot = %s",(notation,))
            
        cases = cur.fetchall()
        

        for case in cases:
            elnot = case[0]
            mod_type = case[1]
            mine_global_clstrs(elnot,mod_type)
            
    ## After mining new clstrs, run residue through the associator
    #  - some residue_ints may now match due to wobulation
    #  - new clstrs are available to attach
    
    

    ##2. form new modes
    if notation == None:
        cur.execute("select elnot, adj_mt(mod_type) mt, count(*) from " +fq_idb+ """intercepts where intercept_id in
                (select intercept_id from residue_ints where r_type='P') group by elnot,mt having count(*) >= %s""",(min_residue_ints,))
    else:
       cur.execute("select elnot, adj_mt(mod_type) mt, count(*) from " +fq_idb+ """intercepts where elnot = %s and intercept_id in
                (select intercept_id from residue_ints where r_type='P') group by elnot,mt having count(*) >= %s""",(notation,min_residue_ints)) 
    cases = cur.fetchall()
    
    for case in cases:
        elnot = case[0]
        mod_type = case[1]
        print(elnot,mod_type)
        
        ### Pre associate residue data (there may be conformal sequences in residue due to wobulation)
        cur.execute("delete from comparison_snapshot")
        cur.execute("commit")
        cur.execute("insert into comparison_snapshot(intercept_id) select a.intercept_id from residue_ints a,"+fq_idb+"""intercepts b
                    where a.intercept_id=b.intercept_id and b.elnot = %s and adj_mt(b.mod_type) = %s""",(elnot,mod_type))
        new_associator(False)
        
        
        ### Extract new sequences from residue
        # Since new clusters already stored just retrieve available pri clusters
        cur.execute("select clstr_id,parm_min,parm_max from global_parm_clstrs where elnot=%s and mod_type=%s and parm='PRI' order by parm_min",
                    (elnot,mod_type))
        pri_clstrs = cur.fetchall()
    
        #retrieve intercept (sequence) residue
        cur.execute("select a.intercept_id,pri_value,pri_number from " +fq_idb+ """intercept_pris a, intercepts b, residue_ints c where
                    a.intercept_id=b.intercept_id and b.intercept_id=c.intercept_id and c.r_type='P' and b.elnot=%s and adj_mt(b.mod_type)=%s""",
                    (elnot,mod_type))
        residue = cur.fetchall()
    
           
    
        # Express intercept sequences in clstr_id lingo
        
        seqs,seq_ints = get_residue_seqs(residue,pri_clstrs,min_seq_cnt)   # seqs: list of [seq_num, [seq], seq_cnt], seq_ints: list of [int_id, seq_num]
                                                                           # seq_num = 1-up index, [seq] is a list of clstr_id (ordered)

        # Group similar sequences into parents & children
        parent_seqs,child_seqs = group_seqs(seqs,mod_type)    # parent_seqs: list of [seq_num,[coded_seq], hc]
                                                              # child_seqs: list of [seq_num (parent), [coded child_seq], child_hc, seq_num (child)]
        #return parent_seqs,child_seqs
    
        # Insert new sequences & sequence_ints
        #  -- Note: child_seqs only explicitly saved if save_children flag = True in VLconfig (otherwise child_ints associated directly to parent)
        #db_bulk_clstr_insert(elnot,mod_type,'PRI',new_pri_clstrs,False)
        new_modes = db_insert_seqs(elnot, mod_type, parent_seqs, child_seqs, seq_ints, False, True)  #last 2 args deletes residue_ints & returns new mode_ids if True
        print(len(new_modes),"new parameter groups found: ",elnot,mod_type)                                                                                            # new_modes: list of mode_id
        print("new_modes",new_modes)
        #### attach clstrs to new modes
        # retrieve relevant clusters
        
        db_clstrs = get_global_clstrs(elnot,mod_type)
        mode_parm_inserts = []
        for mode_id in new_modes:
            # retrieve parm values for associated ints
            # ID global clstr_ids
            # extract local min, local max
            # update local_parm_clstrs

            mode_dat = get_mode_dat(mode_id)
            parm_inserts = get_parm_inserts(db_clstrs,mode_dat)
            mode_parm_inserts.append([mode_id,parm_inserts])
        #print("mode_parm_inserts",mode_parm_inserts)
        
        # insert local parms for all new modes
        db_insert_local_parms(mode_parm_inserts)
                
        # delete residue at end
        cur.execute("delete from residue_ints where r_type='P' and intercept_id in (select intercept_id from sequence_ints)")
        cur.execute("commit")
        
    cur.close()
    db.close()

    
    
# called by new_residue_miner
def get_parm_inserts(clstrs,data):
    # clstrs = [ [rf_clstrs], [pd_clstrs], [sp_slstrs], [ir_clstrs] ]
    # data = [ [rf_vals], [pd_vals], [sp_vals], [ir_vals] ]
    
    i = 0
    parm_inserts = []
    for row in clstrs:   #one row per parm
        dat = data[i]
        if i == 0:
            lparm = 'RF'
        elif i == 1:
            lparm = 'PD'
        elif i == 2:
            lparm = 'SP'
        elif i == 3:
            lparm = 'IR'
        i +=1
        
        for clstr in row:
            clstr_id = clstr[0]
            cmin = clstr[1]
            cmax = clstr[2]
            new_clstr = []
            l_min = -1
            l_max = -1
            for val in dat:
                if val >= cmin and val < cmax:
                    new_clstr = add_unique(clstr_id,new_clstr)
                    if l_min == -1:
                        l_min = val
                    elif val < l_min:
                        l_min = val
                    if l_max == -1:
                        l_max = val
                    elif val > l_max:
                        l_max = val
            if new_clstr != []:
                parm_inserts.append([lparm,clstr_id,l_min,l_max])
    
    #print("parm_inserts",parm_inserts)    
    return parm_inserts
                
                
                    
        

#called by new_residue_miner
def get_global_clstrs(elnot,mod_type):
    db = db_con()
    cur = db.cursor()
    
    cur.execute("select clstr_id,parm,parm_min,parm_max from global_parm_clstrs where elnot = %s and mod_type=%s order by parm,parm_min",(elnot,mod_type))
    db_clstrs = cur.fetchall()
    rf_clstrs = []
    pd_clstrs = []
    sp_clstrs = []
    ir_clstrs = []
    for clstr in db_clstrs:
        parm = clstr[1]
        if parm == 'RF':
            rf_clstrs.append([clstr[0],float(clstr[2]),float(clstr[3])])
        elif parm == 'PD':
            pd_clstrs.append([clstr[0],float(clstr[2]),float(clstr[3])])
        elif parm == 'SP':
            sp_clstrs.append([clstr[0],float(clstr[2]),float(clstr[3])])
        elif parm == 'IR':
            ir_clstrs.append([clstr[0],float(clstr[2]),float(clstr[3])])
    
    cur.close()
    db.close()
    return [rf_clstrs,pd_clstrs,sp_clstrs,ir_clstrs]
    
# called by new_residue_miner
def get_mode_dat(mode_id):
    db = db_con()
    cur = db.cursor()
    cur.execute("select rf,pd,sp,ir,scan_type from " +fq_idb+ """intercepts a, residue_ints b, sequence_ints c, sequences d where 
                a.intercept_id = b.intercept_id and b.intercept_id = c.intercept_id and c.seq_id = d.seq_id and b.r_type = 'P'
                and d.group_id = %s""",(mode_id,))
    dat = cur.fetchall()
    rf_dat = []
    pd_dat = []
    sp_dat = []
    ir_dat = []
    st_dat = []
    for row in dat:
        rf_val = row[0]
        pd_val = row[1]
        sp_val = row[2]
        ir_val = row[3]
        st_val = row[4]
        rf_dat.append(float(rf_val))
        if pd_val != None and pd_val > 0:
            pd_dat.append(float(pd_val))
        if sp_val != None and sp_val > 0:
            sp_dat.append(float(sp_val))
        if ir_val != None and ir_val > 0:
            ir_dat.append(float(ir_val))
        if st_val != None and st_val != '-' and st_val != 'Z':
            st_dat = add_unique(st_val,st_dat)
    st_dat = rip(st_dat,0)
    rf_dat.sort()
    pd_dat.sort()
    sp_dat.sort()
    ir_dat.sort()
    st_dat.sort()
        
    return [rf_dat,pd_dat,sp_dat,ir_dat,st_dat]
    
    cur.close()
    db.close()

# called by new_residue_miner
def mine_global_clstrs(elnot,mod_type):
    
    for parm in ['RF','PRI','PD','SP','IR']:
        print("Mining",parm,"clusters")
        new_clstrs = []
        
        incr,horizon,thresh,p = get_tuning(elnot,parm,'BOTH')
        parm_slope = retrieve_residue(elnot,mod_type,parm,incr)
        if parm_slope == []:
            continue
        
        
        new_clstrs = extract_residue_clstrs(parm_slope,incr,horizon,thresh)
        if new_clstrs == []:
            print("     ...No clstrs extracted from residue")
            continue
        
        # deconflict clusters: if the mined cluster is less than a horizon away from an existing cluster, remove it...
        #                      wobulation is responsible for expanding existing clusters as appropriate
        new_clstrs = deconflict_clusters(elnot,mod_type,parm,new_clstrs,horizon)
        
        ## append bins prior to insert
        clstr_inserts = []
        bin_lims = []
        for clstr in new_clstrs:
            cmin = clstr[0]
            cmax = clstr[1]
            clstr_bins,bmin,bmax = get_bins(parm_slope, cmin, cmax, incr, horizon * k_buf)
            clstr_inserts.append([cmin,cmax,clstr_bins])
            bin_lims.append([bmin,bmax])
        
        
        ## insert new clstrs into db
        db_global_parm_insert(elnot,mod_type,parm,clstr_inserts)
        
        ## purge residue
        k = round(1/incr)
        for row in bin_lims:
            bmin = row[0]
            bmax = row[1]
            remove_residue(elnot, mod_type, parm, k*bmin, k*bmax)
            
        

# called in residue_miner.mine_global_clstrs
def get_bins(parm_slope, clstr_min, clstr_max, incr, buffer):
    k = round(1/incr)
    num_bins = round(((clstr_max + buffer) - (clstr_min - buffer)) * k)
    bin_min = round(clstr_min-buffer,1)
    bin_max = round(clstr_max+buffer,1)
    bins = num_bins * [0]
    for i in range(0,num_bins):
        bin_val = round(bin_min + i*incr,1)
        for row in parm_slope:
            idx = float(row[0])
            if feql(idx,bin_val):
                bins[i] = row[1]
    return bins,bin_min,bin_max
        
            
        
#called in residue_miner.mine_global_clstrs    
def deconflict_clusters(elnot,mod_type,parm,new_clstrs,horizon):
    #print("***** called Deconflict clstrs******")
    db = db_con()
    cur = db.cursor()
    cur.execute("select clstr_id,parm_min,parm_max from global_parm_clstrs where elnot=%s and mod_type=%s and parm=%s order by parm_min",
                (elnot,mod_type,parm))
    dat = cur.fetchall()
    db_clstrs = []
    for row in dat:
        db_clstrs.append([row[0],float(row[1]),float(row[2])])
    new_clstrs = absorb_clstrs(db_clstrs,new_clstrs,horizon)
    cur.close()
    db.close()
    return new_clstrs

# called in residue_miner.deconflict_clusters    
def absorb_clstrs(db_clstrs, new_clstrs, horizon):
    """Global clstr values have already been used to wobulate db clusters, so simply need to 
       remove any new clusters falling within the bins (+/- horizon) of existing clstrs
       db_clstrs: list of [clstr_id, clstr_min, clstr_max]
       new_clstrs: list of [clstr_min, clstr_max]
    """
    #print("******absorb clstrs called")
    for new_clstr in new_clstrs.copy():
        for db_clstr in db_clstrs.copy():
            tst_clstrs = [db_clstr[1:3],new_clstr]
            if len(simple_merge(tst_clstrs,horizon)) == 1:
                new_clstrs.remove(new_clstr)
                break
    return new_clstrs

# called in residue_miner.mine_global_clstrs
def extract_residue_clstrs(parm_slope,incr,horizon,thresh):
    # Since residue is already in the form of the legacy parm_slope, just need to extract clstrs
    
    max_slope = max(rip(parm_slope,1))
    if max_slope < min_clstr_peak/incr:
        print("failed to meet min slope threshold")
        print("parm_slope",parm_slope)
        return []

    out_clstrs = []
    in_clstr = False
    in_lull = False
    avg_slope = 0
    k = 0
    last_valid = 0
    clstr_min = -1
    clstr_count = 0
    clstr_count2 = 0
    
    for val in parm_slope:
        if not in_clstr:
            if val[1] < thresh * max_slope:
                continue
            else:
                in_clstr = True
                clstr_min = float(val[0])
                last_valid = clstr_min
                clstr_count = val[1]
                clstr_count2 = 0
                in_lull = False
                avg_slope = 0
                k = 0
        else:  # in a clstr
            if float(val[0]) - last_valid > horizon:   #cao of prior clstr and start a new one if thresh exceeded or continue search
                clstr_max = last_valid + incr
                if clstr_count-clstr_count2 >= min_clstr_volume/incr:
                    out_clstrs.append([trunc2(clstr_min,incr),trunc2(clstr_max,incr)])
                if val[1] >= thresh * max_slope:
                    clstr_min = float(val[0])
                    last_valid = float(val[0])
                    clstr_count += val[1]
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
                    last_valid = float(val[0])
                    clstr_count += val[1]
                    avg_slope = 0
                    continue
                else: #not in a lull yet but lost threshold
                    clstr_count += val[1]
                    clstr_count2 = val[1]
                    avg_slope = val[1]
                    k = 1
                    in_lull = True
                    continue
            else:  # already in a lull
                clstr_count += val[1]
                clstr_count2 += val[1]
                k = round((float(val[0])-last_valid)/incr)
                avg_slope = ( k*avg_slope + val[1] ) / (k+1)
                if avg_slope >= thresh * max_slope:  
                    clstr_count2 = 0
                    last_valid = float(val[0])
                    in_lull = False
                    avg_slope = 0
                    k = 0
                    
    if in_clstr == True and clstr_count-clstr_count2 >= min_clstr_volume/incr:
        clstr_max = last_valid + incr
        out_clstrs.append([trunc2(clstr_min,incr),trunc2(clstr_max,incr)])

    return out_clstrs

# called in residue_miner.mine_global_clstrs    
def retrieve_residue(elnot,mod_type,parm, incr = 0.1):
    # Retrieves residue in the form of paem_slope
    # parm_residue of form: elnot,mod_type, parm, idx
    
    db = db_con()
    cur = db.cursor()

    k = round(1/incr)
    cur.execute("""select idx* %s, idx_cnt * %s from parm_residue where elnot = %s and mod_type = %s and parm = %s
                    order by idx""",(incr,k,elnot,mod_type,parm))    
    
    residue = cur.fetchall()
    if sum(rip(residue,1)) * incr < min_residue_ints:
        print("insufficient residue available:",elnot,parm,mod_type)
        return []
    else:
        residue.sort()
        return residue
    
    cur.close()
    db.close()
    


# called in residue_miner.mine_global_clstrs
def remove_residue(elnot,mod_type,parm,min_val,max_val):
    # min_val/max_val represent ?  wait till we need it to decide
    #
    db = db_con()
    cur = db.cursor()
    cur.execute("delete from parm_residue where elnot = %s and mod_type = %s and parm = %s and idx between %s and %s",
                (elnot,mod_type,parm,min_val,max_val))
    cur.execute("commit")
    cur.close()
    db.close()

# called in residue_miner.mine_global_clstrs
def db_global_parm_insert(elnot,mod_type,parm,parm_inserts):
    
    db = db_con()
    cur = db.cursor()
    
    cur.execute("select max(clstr_id) from global_parm_clstrs")
    next_id = cur.fetchone()[0]
    if next_id == None:
        next_id = 1
    else:
        next_id = next_id + 1
    for clstr in parm_inserts:
        pmin = clstr[0]
        pmax = clstr[1]
        bins = clstr[2]
        cur.execute("""insert into global_parm_clstrs(clstr_id,elnot,mod_type,parm,parm_min,parm_max,bin_counts)
                       values(%s,%s,%s,%s,%s,%s,%s)""",(next_id,elnot,mod_type,parm,pmin,pmax,bins))
        next_id +=1
    
    cur.execute("commit")
    cur.close()
    db.close()


    
    
