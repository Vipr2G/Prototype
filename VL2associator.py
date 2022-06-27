#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 12:28:35 2021

@author: russelmiller
"""

from VL2db_utils import db_con, get_tuning, log_event, get_db_seqs, add_residue
from VL2dm_utils import sim_seq, encode_int_pris
from VL2config import fq_idb,stage_geo, max_db_trans,max_db_list_len
from VL2shared_utils import rip, trunc2,simple_merge,feql,add_unique
from VL2bin_utils import update_bin_cnts, bin_resize, calc_bin_parms
import operator
import numpy as np

def new_associator(new_data = True):
    # Assumptions:
    #   - data is staged in comparison_queue to be associated
    #   - an environment exists for intercepts to be associated to
    # Function:
    #   - Move data to comparison_snapshot
    #   - Wobulate global clusters (all parms)
    #   - associate ints & new clstrs to existing parameter groups
    #
    db = db_con()
    cur = db.cursor()
    
    ##################################################
    #### Verify there is an ENV to associate to
    
    cur.execute("select count(1) from sequences")
    env_cnt = cur.fetchone()[0]
    if env_cnt == 0:
        event_str = "INFO: Associator called when no ENV existed"
        print(event_str)
        log_event(event_str)
        return
    
    cur.execute("select count(1) from comparison_queue")
    env_cnt = cur.fetchone()[0]
    if env_cnt == 0 and new_data:
        event_str = "INFO: Associator called with empty comparison_queue"
        print(event_str)
        log_event(event_str)
        return

    
    ###################################################
    #### Initialize
    if new_data:      # residue_miner pre-populates comparison_snapshot when pre-processing residue
        cur.execute ("delete from comparison_snapshot")
        cur.execute ("insert into comparison_snapshot(intercept_id) select intercept_id from comparison_queue")
        

        cur.execute("delete from comparison_queue where intercept_id in (select intercept_id from comparison_snapshot)")
        
        
        if stage_geo == True:
            cur.execute("delete from comparison_geo_snapshot")
            cur.execute("""insert into comparison_geo_snapshot(intercept_id) select b.intercept_id from """ + fq_idb + """intercepts a,
                        comparison_snapshot b where a.intercept_id = b.intercept_id and substr(a.elnot,1,1) not in ('A','C','K','M','O')""")
        cur.execute("commit")
        
    
    cur.execute("""select distinct elnot,adj_mt(mod_type) mt from """ + fq_idb + """intercepts a,comparison_snapshot b where 
                    a.intercept_id = b.intercept_id order by elnot,mt""")
   
        
    cases = cur.fetchall()
    
    if len(cases) == 0:
        print("No new data to associate...exiting")
        return
       
    ####################################################
    #### Begin association
    
    for case in cases:
        elnot = case[0]
        mod_type = case[1]
        
        ### wobulate_global clusters:
        #    - wobulates and updates global_parm_clusters
        #    - merges clusters when apropriate and performs housekeeping
        #    - updates parm residue with unmatched data
        #    - commits all changes to DB (subsequent association based exclusively on DB content)
        if new_data:
            wobulate_global_clstrs(elnot,mod_type)
        
        ### Associate ints to existing sequences
        #     Returns
        #     - res_ints: list of int_id
        #     - seq_updates: sorted list of [seq_id,int_id]
        seq_updates,res_ints = associate_ints(elnot,mod_type)     
                                                                  

        ### Update database with the sequence association results:
        ####  - unassociated ints --> residue
        ####  - update sequence_ints

        # insert residue
        if new_data:
            cnt = 0
            for r_int_id in res_ints:
                cur.execute("insert into residue_ints(intercept_id,r_type) values(%s,'P')",(r_int_id,))
                cnt +=1
                if cnt > max_db_trans:
                    cur.execute("commit")
                    cnt = 0
            cur.execute("commit")
            
        # insert seq_ints
        if seq_updates == []:
            continue
        
        cnt = 0
        for seq_update in seq_updates:
            seq_id = seq_update[0]
            int_id = seq_update[1]
            cur.execute("insert into sequence_ints(seq_id,intercept_id) values(%s,%s)",(seq_id,int_id))
            cnt +=1
            if cnt > max_db_trans:
                cur.execute("commit")
                cnt = 0
                
            if not new_data:
                cur.execute("delete from residue_ints where intercept_id = %s and r_type = 'P'",(int_id,))
 
        cur.execute("commit")

        ### Attach/update local_clstrs to updated modes
        # - map associated ints to global clstrs
        # - summarize by mode/by parm each represented clstr and the max/min for all data associated
        #   - if clstr already exists on the mode, update local_min/max if appropriate
        #   - if clstr does not exist on the mode, add it along with the local_min/max
        # - find min/max for each clstr group
        
        # retrieve global clstrs
        cur.execute("""select clstr_id,parm,parm_min,parm_max from global_parm_clstrs where elnot = %s 
                    and mod_type = %s order by parm,parm_min""",(elnot,mod_type))
        global_clstrs = cur.fetchall()
        
        update_local_parms(global_clstrs, seq_updates)
        
    cur.close()
    db.close()



def associate_ints(elnot,mod_type):
    # called by new_associator
    # Assumptions:
    #    - data is already staged in comparison_snapshot
    #    - global clusters have already been wobulated
    # Function:
    #    - Associate staged ints to modes by matching PRI sequence
    #    - Attach any new secondary parameter clusters to existing modes based on associated ints
    #    - Stage unassociated ints in residue_ints (r_type = P)
    db = db_con()
    cur = db.cursor()
    
    #### Retrieve relevant PRI clstrs from the IPE
    cur.execute("select clstr_id,parm_min,parm_max,bin_counts from global_parm_clstrs where parm = 'PRI' and elnot = %s and mod_type = %s order by parm_min",
                (elnot,mod_type))
    pri_clstrs = cur.fetchall()      # Tuple of [clstr_id, pri_min, pri_max, [bins]]
    
    #### Retrieve relevant PRI data
    cur.execute("""select a.intercept_id, a.pri_value, a.pri_number from """ + fq_idb + """intercept_pris a, """ + fq_idb + """intercepts b, comparison_snapshot c where 
                    a.intercept_id = b.intercept_id and a.intercept_id=c.intercept_id and b.elnot = %s and adj_mt(b.mod_type) = %s
                    order by a.intercept_id, pri_number""",(elnot,mod_type))            
    int_pri_dat = cur.fetchall()      # Tuple of [int_id, pri_val, pri_num]
        

    # paste line ----------------------
    ## Express intercept sequences in terms of clstrs
    int_seqs = encode_int_pris(int_pri_dat,pri_clstrs)

    if int_seqs == []:
        res_ints = rip(int_pri_dat,0)
        return [],res_ints
    
    # int_seqs contains only valid seqs and list elements are of the form [int_id, clstr_id1,....,clstr_id_n]  
        
    
    ## Stage ints that couldn't form sequences for residue
    int_seq_ids = rip(int_seqs,0)
    res_ints = []
    prev_id = 0
    for ints in int_pri_dat:
        iid = ints[0]
        if iid != prev_id:
            if iid not in int_seq_ids:
                res_ints.append(iid)
        prev_id = iid
        
                
    # de-dup int seqs into form [[<list of int_ids>], [<list of element_ids>]]  #important to not create duplicate children
    # dedup fctn modified in VL2..."dup" --> exact match req'd, no sim_seq
    unq_int_seqs = dedup_int_seqs(int_seqs,mod_type)
    
        
    #### If valid int_seqs found:
    if len(unq_int_seqs) != 0:
        # retrieve all parent sequences for (elnot,mod_type)
        env_seqs = get_db_seqs(elnot,mod_type,'T')   #retrieve existing parent seqs for elnot, mod_type 
        # list of elements of the form [seq_id, [elem1...elem_n]]
            
        seq_updates = []
        for unq_int_seq in unq_int_seqs:
            res_insert = True
            int_seq = unq_int_seq[1]
            int_ids = unq_int_seq[0]
            for parents in env_seqs:
                parent_seq = parents[1]
                parent_id = parents[0]
                if sim_seq(int_seq,parent_seq,mod_type):
                    for int_id in int_ids:
                        seq_updates.append([parent_id,int_id])
                        #seq_updates = add_unique([parent_id,int_id],seq_updates)
                    res_insert = False
                    # no break to allow ints to match multiple parent seqs
            if res_insert:    # captures additional residue that formed valid seq but didn't match anything in IPE
                for int_id in int_ids:
                    res_ints.append(int_id)
                    
        seq_updates.sort()
    
    cur.close()
    db.close()
    return seq_updates,res_ints

    
    
def dedup_int_seqs(int_seqs,mod_type):
    # called by associator.associate_ints
    # modified for VL2...unique as reported, no sim_seq, etc...that's all dealt with on association
    unique_seqs = []
    int_unq_seqs = []
    if len(int_seqs) == 0:
        return []
    for row in int_seqs:
        int_id = row[0]
        seq = row[1:]
        unique_seqs = add_unique(seq, unique_seqs)
    unique_seqs = rip(unique_seqs,0)
    
    # remove "similar" sequences that are the same length (they are the same seq reported differently)
    #tst_seqs = unique_seqs.copy()
    #for i in range(0,len(unique_seqs)):
    #    sq1 = tst_seqs[i]
    #    for j in range(i+1,len(unique_seqs)):
    #        sq2 = tst_seqs[j]
    #        if sim_seq(sq1,sq2,mod_type) and len(sq1) == len(sq2):
    #            try:
    #                unique_seqs.remove(sq2)
    #            except:
    #                None
                
    for seq in unique_seqs:
        int_list = []
        for row in int_seqs:
            int_id = row[0]
            int_seq = row[1:]
            if int_seq == seq: #or (sim_seq(int_seq,seq,mod_type) and len(int_seq) == len(seq)):
                int_list.append(int_id)
        int_unq_seqs.append([int_list,seq])
    return int_unq_seqs    
    
    
    
def wobulate_global_clstrs(elnot,mod_type):
    # called by new_associator
    db = db_con()
    cur = db.cursor()
    
    for parm in ['PRI','RF','PD','SP','IR']:
        # 1. retrieve new data from database
        if parm == 'PRI':
            print("wobulating PRI...")
            cur.execute("""select pri_value from """ + fq_idb + """intercept_pris a, """ + fq_idb + """intercepts b, comparison_snapshot c where 
                    a.intercept_id = b.intercept_id and a.intercept_id=c.intercept_id and b.elnot = %s and adj_mt(b.mod_type) = %s
                    order by pri_value""",(elnot,mod_type))
        elif parm == 'RF':
            print("wobulating RF....")
            cur.execute("""select rf from """ + fq_idb + """intercepts a, comparison_snapshot b where 
                    a.intercept_id = b.intercept_id and a.elnot = %s and adj_mt(a.mod_type) = %s
                    order by rf""",(elnot,mod_type))
        elif parm == 'PD':
            print("wobulating PD...")
            cur.execute("""select pd from """ + fq_idb + """intercepts a, comparison_snapshot b where 
                    a.intercept_id = b.intercept_id and a.elnot = %s and adj_mt(a.mod_type) = %s and pd > 0
                    order by pd""",(elnot,mod_type))
        elif parm == 'SP':
            print("wopulating SP...")
            cur.execute("""select sp from """ + fq_idb + """intercepts a, comparison_snapshot b where 
                    a.intercept_id = b.intercept_id and a.elnot = %s and adj_mt(a.mod_type) = %s and sp > 0
                    order by sp""",(elnot,mod_type))
        elif parm == 'IR':
            print("wobulating IR...")
            cur.execute("""select ir from """ + fq_idb + """intercepts a, comparison_snapshot b where 
                    a.intercept_id = b.intercept_id and a.elnot = %s and adj_mt(a.mod_type) = %s and ir > 0
                    order by ir""",(elnot,mod_type))
        parm_dat = rip(cur.fetchall(),0)
        
        if parm_dat == []:
            print("    ...NO",parm,"data to wobulate")
            continue
        
        # 2. retrieve relevant db clusters
        cur.execute("""select clstr_id,parm_min,parm_max,bin_counts from global_parm_clstrs where elnot = %s and mod_type = %s 
                    and parm = %s order by parm_min""",(elnot,mod_type,parm))
        db_clstrs = cur.fetchall()
        if db_clstrs ==[]:
            print("    ...no db clstrs exist")
    
            #continue
        
        
        # 3. Wobulate the clusters
        incr,horizon,p = get_tuning(elnot,parm)
        clstr_updates, clstr_deletes,residue = wobulate_all(elnot,parm,parm_dat,db_clstrs,horizon,incr,p)
        
        # 4. Update the db & perform housekeeping
        install_updates2(elnot, parm, mod_type, clstr_updates, clstr_deletes,residue,incr)
            
    
    cur.close()
    db.close()


def wobulate_all(elnot,parm,w_data,w_clstrs,horizon,incr,p):
    # called in associator.wobulate_global_clstrs
    # elnot/parm inputs for feedback only
    wdata = w_data.copy()
    wclstrs = w_clstrs.copy()
    
    if len(wdata) == 0:  # Note: if there are no clstrs need to continue so that parm_residue gets updated
        return [],[],[]
    
    c_updates = []
    clstr_inserts = []
    clstr_deletes = []
    residue = w_data.copy()
    for clstr in wclstrs:
         clstr_id = clstr[0]
         clstr_min = float(clstr[1])
         clstr_max = float(clstr[2])
         clstr_bins = clstr[3]
         bin_min = trunc2(clstr_min-horizon,incr)
         bin_max = trunc2(clstr_max+horizon,incr)
         new_bins = clstr_bins.copy()
         
         # Update bin counts with new data
         for val in wdata:
             prev_bins = new_bins
             new_bins = update_bin_cnts(val,new_bins,clstr_min,clstr_max,horizon,incr)
             if new_bins != prev_bins:
                 try:
                     residue.remove(val)
                 except:
                     None
         
         # Recalculate clstr_min, clstr_max based on the new bin counts
         #  -- if there is no data corresponding to the current clstr the bins haven't changed and no update is required
         if new_bins != clstr_bins:
             new_clstr_min,new_clstr_max = calc_bin_parms(new_bins,bin_min,incr,p)    # trunc2 floats
             new_bin_min = trunc2(new_clstr_min - horizon,incr)
             new_bin_max = trunc2(new_clstr_max + horizon,incr)
             new_bins = bin_resize(new_bins,incr,bin_min,new_bin_min,bin_max,new_bin_max)
             c_updates.append([clstr_id,new_clstr_min,new_clstr_max,new_bins])   # list of wobulated clusters
    
    #### Validation Section: Ensure updates do not result in overlapping clusters
    test_clstrs = []
    for db_clstr in wclstrs:   # form a complete clstr list with either original boundaries or the updated boundaries
        db_id = db_clstr[0]
        old_min = float(db_clstr[1])
        old_max = float(db_clstr[2])
        update_flag = False
        for mod_clstrs in c_updates:
            mod_id = mod_clstrs[0]
            new_min = float(mod_clstrs[1])
            new_max = float(mod_clstrs[2])
            if db_id == mod_id:
                test_clstrs.append([new_min,new_max])
                update_flag = True
                break
        if update_flag == False:
            test_clstrs.append([old_min, old_max])


    len1 = len(test_clstrs)   #210928 Changed from wclstrs to test_clstrs...shouldn't matter but easier to follow logic
    len2 = len(simple_merge(test_clstrs,horizon))
    if len1 != len2:
        #c_updates, clstr_inserts, clstr_deletes = reconcile_updates(wclstrs, c_updates, incr, horizon, p)
        c_updates,clstr_deletes = reconcile2(wclstrs,c_updates,incr,horizon,p)
        event_str = "INFO (wobulate_all): " + str(len(c_updates)) + " " + parm + " clstrs updated, " + str(len(clstr_inserts)) 
        event_str = event_str + " inserted, and " + str(len(clstr_deletes)) + " deleted: " + elnot + ", "
        log_event(event_str)
        
    return c_updates,clstr_deletes,residue




def reconcile2(db_clstrs,clstr_updates,incr,horizon,p):
    # called in associator.wobulate_all
    # 1. identify which clstrs require merging
    #    a. form post-update clstr list: [clstr_id, clstr_min, clstr_max]

    updated_clstrs = []     # list of [clstr_id,clstr_min, clstr_max] with either original or updated limits
    for clstr in db_clstrs:
        match = False
        clstr_id = clstr[0]
        old_min = clstr[1]
        old_max = clstr[2]
        old_bins = clstr[3]
        for delta in clstr_updates:
            delta_id = delta[0]
            new_min = delta[1]
            new_max = delta[2]
            new_bins = delta[3]
            if delta_id == clstr_id:
                updated_clstrs.append([clstr_id,new_min,new_max,new_bins])
                match = True
                break
        if match == False:
            updated_clstrs.append([clstr_id,old_min,old_max,old_bins])
    
    # 2. determine which clstrs need to be merged
    sorted_clstrs = sorted(updated_clstrs, key=operator.itemgetter(1))
    
    base_idx = 0
    more_clstrs = True
    test_clstrs = []
    merge_groups = []
    while more_clstrs:
        if base_idx >= len(sorted_clstrs)-1:
            break
        test_clstrs = [sorted_clstrs[base_idx][:-1]]
        merge_list = [test_clstrs[0][0]]
        merge = True
        j = 1
        while merge:
            if base_idx + j > len(sorted_clstrs)-1:
                more_clstrs = False
                break
            next_clstr = sorted_clstrs[base_idx + j]
            test_clstrs.append(next_clstr[:-1])
            if len(simple_merge(test_clstrs,horizon)) == 1:
                merge_list.append(next_clstr[0])
                j +=1
                continue
            base_idx = base_idx + j
            merge = False
        if len(merge_list) > 1:
            merge_groups.append(merge_list)
    #return merge_groups           
        
    # 3. Merge the clstrs in each merge group
    #    - use first_id for the new clstr_id & pass merge groups to update new_id for absorbed clstrs in db
        

    
    merged_clstrs = []
    for row in merge_groups:
        merge_dat = []
        mid = row[0]
        for cid in row:
            for clstr in updated_clstrs:
                if clstr[0] == cid:
                    merge_dat.append(clstr)
                    break
        new_min,new_max,new_bins = merge_clstrs(merge_dat,incr,horizon,p)
        merged_clstrs.append([mid,new_min,new_max,new_bins])
    # At this point:
    #  merged clstrs have original clstr_ids & data to be updated
    #    -recombine with unmerged updates,,,use "if clstr_id not in obe_ids" ...it works!!
    #  merge_groups hold data for housekeeping (clstrs 2,... are all absorbed into the first value on each row)
    merge_list = []
    for row in merge_groups:
        for cid in row:
            merge_list.append(cid)
    for clstr in updated_clstrs:
        cid = clstr[0]
        if cid not in merge_list:
            merged_clstrs.append(clstr)
    return merged_clstrs,merge_groups


def merge_clstrs(clstrs,incr, horizon,p):
    # called in associator.reconcile2
    """ 
    clstrs list of [clstr_id,clstr_min,clstr_max,bin_counts]
    for PD,SP,IR mode_id required
    return delete_clstrs (list of clstr_ids), and a single insert_clstr: ([clstr_min,clstr_max,clstr_bins])
    
    after merging the bins, the updated bin will need to be stabilized
    """
    clstr_mins = rip(clstrs,1)
    clstr_maxs = rip(clstrs,2)
    new_clstr_min = float(min(clstr_mins))
    new_bin_min = trunc2(new_clstr_min - horizon,incr)
    new_clstr_max = float(max(clstr_maxs))
    new_bin_max = trunc2(new_clstr_max + horizon,incr)
    num_new_bins = round((new_bin_max-new_bin_min)/incr)
    
    
    resized_bins = []
    
    for clstr in clstrs:
        clstr_min = float(clstr[1])
        clstr_max = float(clstr[2])
        clstr_bins = clstr[3]
        bin_min = trunc2(clstr_min-horizon,incr)    #ToDo verify correct convention of bin_max
        bin_max = trunc2(clstr_max + horizon,incr)
        new_clstr_bins = bin_resize(clstr_bins,incr,bin_min,new_bin_min,bin_max,new_bin_max)
        
        # Shouldn't have to do this, but on rare occasion bins come in with incorrect sizes (off by 1-3)...need to find root cause
        if len(new_clstr_bins) != num_new_bins:
            event_str = "WARNING: incorrect bin lengths had to be 'fixed' in env_maint.merge_clstrs"
            log_event(event_str)
            b_error = num_new_bins - len(new_clstr_bins)
            if b_error > 0 :   # need to add bins
                if b_error == 1:
                    lpad = 1
                    rpad = 0
                elif b_error % 2 == 0:
                    lpad = round(b_error/2)
                    rpad = lpad
                else:
                    lpad = int(b_error/2)
                    rpad = b_error - lpad
                l_add = [0] * lpad
                r_add = [0] * rpad
                new_clstr_bins = l_add + new_clstr_bins + r_add
            else:   # need to remove bins
                if b_error == -1:
                    rpad = b_error
                    lpad = 0
                elif b_error % 2 == 0:
                    lpad = round(b_error/2)
                    rpad = lpad
                else:
                    rpad = round(b_error/2)
                    lpad = b_error - rpad
                
                new_clstr_bins = new_clstr_bins[-lpad:rpad]
                
        
        resized_bins.append(new_clstr_bins)
    
    updated_bins = []
    for i in range(0,len(resized_bins[0])):
        bin_val = new_bin_min + i * incr
        overlap_cnt = 0
        for j in range(0,len(clstr_mins)):
            if float(clstr_mins[j])-horizon <= bin_val and float(clstr_maxs[j]) + horizon > bin_val:   
                overlap_cnt +=1
        bin_cnt = 0
        for j in range(0,len(clstr_mins)):
            bin_cnt += resized_bins[j][i]
        if overlap_cnt != 0:
            bin_cnt = round(bin_cnt/overlap_cnt)
        else:
            bin_cnt = 0
        updated_bins.append(bin_cnt)

    #return new_clstr_min,new_clstr_max,updated_bins            

    # add for VL2 not previously in env_maint
    # Recharacterize the merged clstrs based on the ASP algorithm
    
    new_cmin,new_cmax,new_bins = stabilize_binned_clstrs(new_clstr_min,new_clstr_max,updated_bins,new_bin_min,new_bin_max,incr,horizon,p)    
    return new_cmin,new_cmax,new_bins



def stabilize_binned_clstrs(cmin,cmax,bins,bin_min,bin_max,incr,horizon,p):
    # called in associator.merge_clstrs
    cmin2,cmax2 = calc_bin_parms(bins,bin_min,incr,p)
    if feql(cmin,cmin2) and feql(cmax,cmax2):
        #unstable = False
        return cmin,cmax,bins
    else:
        unstable = True
        new_bins = bins
        new_bin_min = bin_min
        prev_cmin = cmin
        prev_cmax = cmax
    i = 0
    while unstable:
        i +=1
        new_cmin,new_cmax = calc_bin_parms(new_bins,new_bin_min,incr,p)
        if feql(new_cmin,prev_cmin) and feql(new_cmax,prev_cmax):
            unstable = False
        else:
            new_bin_min = trunc2(new_cmin - horizon, incr)
            new_bin_max = trunc2(new_cmax + horizon,incr)
            new_bins = bin_resize(bins,incr,bin_min,new_bin_min,bin_max,new_bin_max)
            prev_cmin = new_cmin
            prev_cmax = new_cmax
        if i>20:
            event_str = "****WARNUNG: Bins failed to stabilize (VL2associator.stabilize_binned_clstrs)"
            print(event_str)
            log_event(event_str)
            break
    event_str = "bins  SUCCESSFULLY stabilited (stabilize_binned_clstrs"
    log_event(event_str)
    return new_cmin,new_cmax,new_bins


def install_updates2(elnot,parm,mod_type,clstr_updates,clstr_deletes,residue,incr):
    # called in associator.wobulate_global_clstrs
    db = db_con()
    cur = db.cursor()
    
    # update changed clstrs       
    for clstr in clstr_updates:
        cid = clstr[0]
        cmin = clstr[1]
        cmax = clstr[2]
        cbins = clstr[3]
        cur.execute("""update global_parm_clstrs set parm_min = %s,parm_max = %s,bin_counts = %s where clstr_id = %s and elnot = %s and
                    mod_type = %s and parm = %s""",(cmin,cmax,cbins,cid,elnot,mod_type,parm))
    cur.execute("commit")
    print(len(clstr_updates),"global",parm,"clstrs updated")
    
    # delete clstrs absorbed in a merge
    for row in clstr_deletes:
        new_id = row[0]
        cnt = 0
        for clstr in row[1:]:
            cnt +=1
            obe_id = clstr
            
            #Since global_clstr_ids are unique, the below corrections *should* be ok
            
            #cur.execute("""update local_parm_clstrs a, parameter_groups b set a.clstr_id = %s where a.group_id=b.group_id and
            #            a.clstr_id = %s and a.parm = %s and b.elnot = %s and b.mod_type = %s""",(new_id,obe_id,parm,elnot,mod_type))
            

            cur.execute("""update local_parm_clstrs set clstr_id = %s where clstr_id = %s """,(new_id,obe_id))
                        
                        
            #cur.execute("""delete from global_parm_clstrs where clstr_id = %s and elnot = %s and mod_type = %s and parm = %s""",
            #            (obe_id,elnot,mod_type,parm))
            cur.execute("delete from global_parm_clstrs where clstr_id = %s",(obe_id,))
            
            if parm == 'PRI':
                #cur.execute("""update sequence_elements a,sequences b, parameter_groups c set a.clstr_id = %s where
                #           a.clstr_id = %s and a.seq_id=b.seq_id and b.group_id = c.group_id and c.elnot = %s and b.mod_type = %s""",
                #            (new_id,obe_id,elnot,mod_type))
                cur.execute("select distinct seq_id from sequence_elements where clstr_id = %s",(obe_id,))
                seq_str = str(cur.fetchall())
            
                event_str = "PRI cldtr_id " + str(obe_id) + " replaced by " +str(new_id) +" seqs " +seq_str
                print(event_str)
                log_event(event_str)
                cur.execute("""update sequence_elements set clstr_id = %s where clstr_id = %s""",(new_id,obe_id))
                cur.execute("commit")
            # Note, this may result in duplicate sequences...issue?
    cur.execute("commit")
    
    # update residue
    add_residue(elnot,mod_type,parm,residue,incr)
    cur.close()
    db.close()
        
def update_local_parms(global_clstrs, seq_updates):
    # called by new_associator
    db = db_con()
    cur = db.cursor()
    
    clstr_inserts = []
    clstr_updates = []
        
    # Get list of updated sequences:
    unq_seq_ids = list(np.unique(rip(seq_updates,0)).astype(int))
    
    # Get list of associated int_ids:
    associated_ints = list(np.unique(rip(seq_updates,1)))
    ints_len = len(associated_ints)
    
    for parm in ['RF','PD','SP','IR']:
        
        # collect global clstrs
        parm_clstrs = list(i for i in global_clstrs if i[1] == parm)
        
        # collect intercept_values
        if len(associated_ints) < max_db_list_len:
            if len(associated_ints) > 1:
                int_str = str(tuple(associated_ints))
            else:
                int_str = "(" + str(associated_ints[0]) + ")"
            cur.execute("select intercept_id, " +parm+ " from " +fq_idb+ "intercepts where " +parm+ """ > 0 
                        and intercept_id in """ +int_str)
            parm_dat = cur.fetchall()
        else:
            idx = 0
            parm_dat = []
            while idx < ints_len:
                pad = min(max_db_list_len,ints_len-idx)
                int_list = tuple(associated_ints[idx:idx+pad])
                if len(int_list) > 1:
                    int_str = str(int_list)
                else:
                    int_str = "(" + str(int_list[0]) + ")"
                cur.execute("select intercept_id," +parm+ " from " +fq_idb+ "intercepts where " +parm+ """ > 0 
                            and intercept_id in """ +int_str)
                parm_dat = parm_dat + cur.fetchall()
                idx = idx + pad
        # Express intercept values in terms of global clstrs
        int_parm_clstrs = []
        for row in parm_dat:
            int_id = row[0]
            val = row[1]
            for clstr in parm_clstrs:
                if val >= clstr[2] and val < clstr[3]:
                    cid = clstr[0]
                    int_parm_clstrs.append([int_id,parm,cid,val])
                    break
        
        # For each updated mode, summarize the represented parm_clstrs
        for u_id in unq_seq_ids:
            associated_mode_clstrs = []
            seq_id = int(u_id)
            mode_ints = list((i[1] for i in seq_updates if i[0]==seq_id))
            mode_clstrs = list((i for i in int_parm_clstrs if i[0] in mode_ints))
            unq_mode_clstrs = list(np.unique(rip(mode_clstrs,2)))

            for clstr in unq_mode_clstrs:
                clstr_vals = list((i[3] for i in mode_clstrs if i[2]==clstr))
                lmin = min(clstr_vals)
                lmax = max(clstr_vals)
                associated_mode_clstrs.append([seq_id,clstr,lmin,lmax])
            
            # Retrieve stored local clstrs for comparison
            cur.execute("""select clstr_id, a.group_id, parm_min,parm_max from local_parm_clstrs a, sequences b where
                        a.group_id = b.group_id and b.seq_id = %s and a.parm='""" +parm+"'",(seq_id,))
            stored_mode_clstrs = cur.fetchall()
            stored_mode_cids = rip(stored_mode_clstrs,0)
            for clstr in associated_mode_clstrs:
                if clstr[1] not in stored_mode_cids:
                    seq_id = clstr[0]
                    cur.execute("select group_id from sequences where seq_id = %s",(seq_id,))
                    group_id = cur.fetchone()[0]
                    #group_id = stored_mode_clstrs[0][1]
                    clstr_inserts.append([clstr[1],group_id,parm,clstr[2],clstr[3]])
                else:
                    associated_cid = clstr[1]
                    stored_clstr = list((i for i in stored_mode_clstrs if i[0]==associated_cid))
                    if len(stored_clstr) > 1:
                        print("WARNING:  Trouble in paradise")
                        print("stored_clstr:",stored_clstr)
                    stored_clstr = stored_clstr[0]
                    stored_min = stored_clstr[2]
                    stored_max = stored_clstr[3]
                    associated_min = clstr[2]
                    associated_max = clstr[3]
                    new_min = min(stored_min,associated_min)
                    new_max = max(stored_max, associated_max)
                    if new_min < stored_min or new_max > stored_max:
                        group_id = stored_clstr[1]
                        clstr_updates.append([associated_cid,group_id,parm,new_min,new_max])
    # Update the database
    cnt = 0
    for row in clstr_updates:
        clstr_id = int(row[0])
        group_id = row[1]
        parm = row[2]
        lmin = row[3]
        lmax = row[4]
        cur.execute("""update local_parm_clstrs set parm_min = %s,parm_max = %s where clstr_id = %s and group_id = %s and parm = %s""",
                    (lmin,lmax,clstr_id,group_id,parm))
        cnt +=1
        if cnt > max_db_trans:
            cur.execute("commit")
            cnt = 0
    
    cnt = 0
    for row in clstr_inserts:
        clstr_id = int(row[0])
        group_id = row[1]
        parm = row[2]
        lmin = row[3]
        lmax = row[4]
        cur.execute("insert into local_parm_clstrs(clstr_id,group_id,parm,parm_min,parm_max) values(%s,%s,%s,%s,%s)",
                    (clstr_id,group_id,parm,lmin,lmax))
        cnt +=1
        if cnt > max_db_trans:
            cur.execute("commit")
            cnt = 0
    
    cur.execute("commit")
    cur.close()
    db.close()


