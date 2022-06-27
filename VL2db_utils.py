#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 12:28:43 2021

@author: russelmiller
"""

#from VL2credentials import user_name, user_password, db_host, db_port, db_name
from VL2config import def_rf_horizon, def_pri_horizon, def_pd_horizon, def_sp_horizon, def_ir_horizon
from VL2config import def_rf_thresh,  def_pri_thresh,  def_pd_thresh,  def_sp_thresh,  def_ir_thresh
from VL2config import def_rf_p,       def_pri_p,       def_pd_p,       def_sp_p,       def_ir_p
from VL2config import def_incr, def_roa, def_grid_res
from VL2config import save_children, fq_idb, max_db_list_len
from VL2shared_utils import add_unique, rip
from math import trunc

import psycopg2
import configparser
config = configparser.ConfigParser()

###############################################
## Externally called fctns
###############################################


def db_con():
    #import psycopg2
    config.read('db.config')
    
    user_name = config['viprlite']['user']
    user_password = config['viprlite']['password']
    db_host = config['viprlite']['host']
    db_port = config['viprlite']['port']
    db_name = config['viprlite']['database']

    db = psycopg2.connect(user = user_name,
                          password = user_password,
                          host = db_host,
                          port = db_port,
                          database = db_name)
    return(db)

def delete_env(elnot, PorG = 'PG', delete_clstrs = True):
    """ Deletes preservet environment data for specified elnot
        delete_env(0) is a global delete of the entire environment for all ELNOTs
    """
    local_db = db_con()
    cur = local_db.cursor()

    cur.execute("delete from comparison_queue")
    cur.execute("delete from comparison_snapshot")
    cur.execute("delete from comparison_geo_snapshot")
    
    if elnot == 0:
        cur.execute("delete from residue_ints")
        cur.execute("delete from sequence_ints")
        cur.execute("delete from sequence_elements")
        cur.execute("delete from sequences")
        cur.execute("delete from local_parm_clstrs")
        
        #cur.execute("delete from wob_modes")
        cur.execute("delete from parameter_groups")
        
        #cur.execute("delete from wob_mode_rfs")
        #cur.execute("delete from wob_mode_sts")
        #cur.execute("delete from wob_pri_clstrs")
        #cur.execute("delete from wob_rf_clstrs")
        cur.execute("delete from wob_event_log")
        #cur.execute("delete from rf_residue_ints")
        cur.execute("delete from global_parm_clstrs")
        cur.execute("delete from parm_residue")
        
        if PorG.upper() == 'PG' or PorG.upper() == 'GP' or PorG.upper() == 'G':
            cur.execute("delete from site_geos")
            cur.execute("delete from site_ints")
    
    else:
        cur.execute("""delete from residue_ints where intercept_id not in 
                       (select intercept_id from intercepts)""")
        cur.execute("""delete from sequence_ints where intercept_id not in 
                       (select intercept_id from intercepts)""")
        #cur.execute("""delete from rf_residue_ints where intercept_id not in 
        #               (select intercept_id from intercepts)""")
                       
        cur.execute("""delete from residue_ints where intercept_id in
                    (select intercept_id from """ + fq_idb + """intercepts where elnot = %s)""",(elnot,))
        
        #cur.execute("""delete from rf_residue_ints where intercept_id in
        #            (select intercept_id from """ + fq_idb + """intercepts where elnot = %s)""",(elnot,))
                    
        cur.execute("""delete from sequence_ints where seq_id in 
                    (select seq_id from sequences where elnot = %s)""",(elnot,))
                     
        cur.execute("""delete from sequence_elements where seq_id in 
                   (select seq_id from sequences where elnot = %s)""",(elnot,))
                    
        cur.execute("delete from sequences where elnot = %s",(elnot,))
        
        cur.execute("""delete from local_parm_clstrs where group_id in
                    (select group_id from parameter_groups where elnot = %s)""",(elnot,))
                    
        
        #cur.execute("""delete from wob_mode_rfs where mode_id in
        #            (select group_id from parameter_groups where elnot = %s)""",(elnot,))
                    
        #cur.execute("""delete from wob_mode_sts where mode_id in 
        #            (select group_id from parameter_groups where elnot = %s)""",(elnot,))
        
        #cur.execute("delete from wob_pri_clstrs where elnot = %s",(elnot,))
        
        #cur.execute("delete from wob_rf_clstrs where elnot = %s",(elnot,))
        
        
        cur.execute("delete from parameter_groups where elnot = %s",(elnot,))
        
        if delete_clstrs:
            cur.execute("delete from global_parm_clstrs where elnot = %s",(elnot,))
            #cur.execute("delete from local_parm_clstrs where elnot = %s",(elnot,))
        
            cur.execute("delete from parm_residue where elnot = %s",(elnot,))
        
            
        if PorG.upper() == 'PG' or PorG.upper() == 'GP' or PorG.upper() == 'G':
            
            cur.execute("""delete from site_ints where site_id in 
                        (select site_id from site_geos where elnot = %s)""",(elnot,))
            cur.execute("delete from site_geos where elnot = %s",(elnot,))            
  
    cur.execute("commit")
    cur.close()
    local_db.close()
    
    
    
def get_tuning(elnot,parm,mode = 'WOB'):
    db = db_con()
    cur = db.cursor()
    
    # Define default values

    # def_incr = 0.1 - imported from config
    if parm == 'RF' or parm == 'RFL':   
        def_horizon = def_rf_horizon     # Min separation between consecutive RF clusters (MHz) 
        def_thresh = def_rf_thresh
        def_p = def_rf_p
    elif parm == 'PRI':
        def_horizon = def_pri_horizon
        def_thresh = def_pri_thresh
        def_p = def_pri_p
    elif parm == 'PD':
        def_horizon = def_pd_horizon
        def_thresh = def_pd_thresh
        def_p = def_pd_p
    elif parm == 'SP':
        def_horizon = def_sp_horizon
        def_thresh = def_sp_thresh
        def_p = def_sp_p
    elif parm == 'IR':
        def_horizon = def_ir_horizon
        def_thresh = def_ir_thresh
        def_p = def_ir_p
    else:
        log_str = "***** WARNING: get_tuning called with invalid parm: " + elnot + " - " + parm
        print(log_str)
        log_event(log_str)
        return
        
    cur.execute("select incr,horizon,p,thresh from direct_tune_parms where elnot = %s and parm = %s",(elnot.upper(),parm.upper()))
    parms = cur.fetchone()
    if parms == None:
        if mode == 'WOB':
            return def_incr, def_horizon, def_p
        elif mode == 'HISTO':
            return def_incr, def_horizon, def_thresh
        elif mode == 'BOTH':
            return def_incr, def_horizon, def_thresh, def_p
        else:
            log_str = "***** WARNING: get_tuning called with invalid mode: " + elnot + " - " + parm + " - " + mode
            print(log_str)
            log_event(log_str)
            return
    else:
        if parms[0] == None:
            incr = def_incr
        else:
            incr = float(parms[0])
            
        if parms[1] == None:
            horizon = def_horizon
        else:
            horizon = float(parms[1])
            
        if parms[2] == None:
            p = def_p
        else:
            p = float(parms[2])
        
        if parms[3] == None:
            thresh = def_thresh
        else:
            thresh = float(parms[3])
            
    if mode == 'WOB':
        return incr,horizon,p
    elif mode == 'HISTO':
        return incr,horizon,thresh
    elif mode == 'BOTH':
        return incr,horizon,thresh,p
    else:
        log_str = "***** WARNING: get_tuning called with invalid mode: " + elnot + " - " + parm + " - " + mode
        print(log_str)
        log_event(log_str)


# VL2 calls: baseline_builder & associator
def add_residue(elnot,mod_type,parm,vals,incr = 0.1):
    # incorporate a list of actual values into residue
    # assumes parm_residue of form: elnot,mod_type, parm, idx, idx_cnt
    
    if len(vals) == 0:
        return
    
    db = db_con()
    cur = db.cursor()
    
    k = round(1/incr)
    unq_idxs = []
    for val in vals:
        l_idx = trunc(val*k)
        unq_idxs = add_unique(l_idx,unq_idxs)
        
    trans_count = 0
    for row in unq_idxs:
        l_idx = row[0]
        l_cnt = row[1]
        cur.execute("select idx_cnt from parm_residue where elnot = %s and mod_type = %s and parm = %s and idx = %s",
                    (elnot,mod_type,parm,l_idx))
        db_cnt = cur.fetchone()
        
        
        if db_cnt != None:
            cur.execute("update parm_residue set idx_cnt = %s where elnot = %s and mod_type = %s and parm = %s and idx = %s",
                        (db_cnt[0]+l_cnt,elnot,mod_type,parm,l_idx))
        else:
            cur.execute("insert into parm_residue(elnot,mod_type,parm,idx,idx_cnt) values(%s,%s,%s,%s,%s)",
                        (elnot,mod_type,parm,l_idx,l_cnt))
        trans_count +=1
        if trans_count > 100:
            cur.execute("commit")
            trans_count = 0

    cur.execute("commit")
    cur.close()
    db.close()
        
# VL2 called from residue_miner
def db_insert_local_parms(parm_inserts):
    db = db_con()
    cur = db.cursor()
    for row in parm_inserts:
        mode_id = row[0]
        clstrs = row[1]
        for dat in clstrs:
            parm = dat[0]
            clstr_id = dat[1]
            lmin = dat[2]
            lmax = dat[3]
            cur.execute("insert into local_parm_clstrs(clstr_id,group_id,parm,parm_min,parm_max) values(%s,%s,%s,%s,%s)",
                        (clstr_id,mode_id,parm,lmin,lmax))
        
    cur.execute("commit") 
    cur.close()
    db.close()


    
    
    
    

def db_insert_seqs(elnot,mod_type, parent_seqs, child_seqs, seq_ints, delete_residue = True,return_modes = False, skip_ints = False, src = 'Residue Miner'):
    db = db_con()
    cur = db.cursor()
    seq_int_inserts = []
    new_modes = []
    cur.execute("select max(seq_id) from sequences")
    seq_id = cur.fetchone()[0]
    if seq_id == None:
        seq_id = 1
    else:
        seq_id +=1
    for parent in parent_seqs:
        p_indx = parent[0]
        p_hc = parent[2]
        # ToDo: consider if need to check DB for an association candidate here b4 insert, if return False then do update as below, else associate internally, e.g.,
        #    db_merge_flag = DB_merge_parent(sq1,sq1_hc)....if False, as is otherwise the fctn takes care of association of both parent and children
        for seq_int in seq_ints:
            if seq_int[1]==p_indx:
                #seq_int_inserts.append([seq_int[0],seq_id])
                seq_int_inserts = add_unique([seq_int[0],seq_id],seq_int_inserts)
                
        # Insert new mode for each parent & retrieve mode_id then insert sequence & sequence elements
        mode_id = db_insert_mode(elnot,mod_type,src)
        new_modes.append(mode_id)
        cur.execute("""insert into sequences(seq_id,elnot,mod_type,group_id,heard_count,is_parent) values(%s,%s,%s,%s,%s,'T')""",
                    (seq_id,elnot,mod_type,mode_id,p_hc))
        # insert sequence_elements (parent)
        i = 0
        for elem in parent[1]:
            clstr_id = elem
            i +=1
            cur.execute("insert into sequence_elements(seq_id,clstr_id,position) values(%s,%s,%s)",(seq_id,clstr_id,i))
        
        #
        #
        ## if not save_children:
        #     just add child int_ids to seq_int_inserts for the appropriate parent & skip bulk of below, else do below
        #
        if not save_children:
            for child in child_seqs:
                c_indx = child[0]
                c_int_key = child[3]
                if c_indx == p_indx:
                    for seq_int in seq_ints:
                        if seq_int[1] == c_int_key:
                            #seq_int_inserts.append([seq_int[0],seq_id])
                            seq_int_inserts = add_unique([seq_int[0],seq_id],seq_int_inserts)
            seq_id +=1
        else:
        
            # insert any child sequences & intercept associations
            pid = seq_id
            seq_id +=1
            for child in child_seqs:
                c_indx = child[0]
                c_hc = child[2]
                c_int_key = child[3]
                if c_indx == p_indx:
                    #k = 0
                    for seq_int in seq_ints:
                        if seq_int[1]==c_int_key:  #indx: # and k > 0:  # Note first match is the parent
                            #seq_int_inserts.append([seq_int[0],seq_id])
                            seq_int_inserts = add_unique([seq_int[0],seq_id],seq_int_inserts)
                        #k +=1
                
                    cur.execute("""insert into sequences(seq_id,elnot,mod_type,group_id,heard_count,is_parent,parent_id) values(%s,%s,%s,%s,%s,'F',%s)""",
                        (seq_id,elnot,mod_type,mode_id,c_hc,pid))
                
                    # insert sequence_elements (child)
                    #ToDo
                    i = 0
                    for elem in child[1]:
                        clstr_id = elem
                        i +=1
                        cur.execute("insert into sequence_elements(seq_id,clstr_id,position) values(%s,%s,%s)",(seq_id,clstr_id,i))
                    seq_id +=1
        cur.execute("commit")
        
        #
        # end legacy indent
        #
        #
    # insert sequence_ints
    if not skip_ints:   #merge related houseKeeping in associator needs to skip intercept inserts cuz it will be taken care of elsewhere
        seq_int_inserts = rip(seq_int_inserts,0)
        for ints in seq_int_inserts:
            cur.execute("insert into sequence_ints(intercept_id,seq_id) values(%s,%s)",(ints[0],ints[1]))
            
    # Clean up residue:
    if delete_residue == True:
        cur.execute("delete from residue_ints where r_type='P' and intercept_id in (select intercept_id from sequence_ints)")
    cur.execute("commit")
    cur.close()
    db.close()
    if return_modes == True:
        return new_modes
##### End db_insert_seqs()

def db_insert_mode(elnot,mod_type,source):
    # called by db_insert_seqs
    #called by db_insert_seqs
    db = db_con()
    cur = db.cursor()
    cur.execute("select max(group_id) from parameter_groups")
    next_id = cur.fetchone()[0]
    if next_id == None:
        next_id = 1
    else:
        next_id +=1
    cur.execute("insert into parameter_groups(group_id,elnot,mod_type,source) values(%s,%s,%s,%s)",(next_id,elnot,mod_type,source))
    cur.execute("commit")
    cur.close()
    db.close()
    return next_id

def get_db_seqs(elnot,mod_type,seq_type = '%'):
    # Return format is [ [int,list],[int,list],...], i.e., [ [seq_id, [elem1,elem2,...]], [seq_id2, [elem1,elem2,...]],...]
    # optional 3rd arg = 'T' : retrieve only parent sequences
    #                    'F' : retrieve only child sequences
    #               no value : return all sequences
    
    db = db_con()
    cur = db.cursor()
    out_seqs = []
    cur.execute("""select min(a.seq_id) from sequence_elements a, sequences b where a.seq_id = b.seq_id and
                b.elnot = %s and b.mod_type = %s and b.is_parent like %s""",(elnot,mod_type,seq_type))
    prev_seq_id = cur.fetchone()[0]
    cur.execute("""select a.seq_id,clstr_id,position from sequence_elements a, sequences b where a.seq_id = b.seq_id and
                b.elnot = %s and b.mod_type = %s and b.is_parent like %s order by a.seq_id,position""",(elnot,mod_type,seq_type))
    cur_seq = []
    for row in cur:
        seq_id = row[0]
        pval = row[1]
        if seq_id == prev_seq_id:
            cur_seq.append(pval)
        else:
            out_seqs.append([prev_seq_id,cur_seq])
            prev_seq_id = seq_id
            cur_seq = [pval]
    # capture last sequence
    out_seqs.append([prev_seq_id,cur_seq])
    return out_seqs
    

def tune_geo(elnot):
    db = db_con()
    cur = db.cursor()
    cur.execute("select count(1) from geo_tuning_parms where elnot = %s",(elnot,))
    parm_rows = cur.fetchone()[0]
    if parm_rows>1:
        event_str = "WARNING: multiple entries for " +elnot+ " in geo_tuning_parms table"
        log_event(event_str)
    elif parm_rows == 1:
        cur.execute("select roa,grid_res from geo_tuning_parms where elnot = %s",(elnot,))
        parms = cur.fetchone()
        roa = parms[0]
        if roa == None:
            roa = def_roa
        else:
            roa = float(roa)
        grid_res = parms[1]
        if grid_res == None:
            grid_res = def_grid_res
        else:
            grid_res = float(grid_res)
    else:
        roa = def_roa
        grid_res = def_grid_res
    cur.close()
    db.close()
    return grid_res, roa
    
    
def db_geo_insert(elnot,sites,site_ints):
    if len(sites) == 0 or len(site_ints) == 0:
        return
    db = db_con()
    cur = db.cursor()
    cur.execute("select max(site_id) from site_geos")
    max_id =cur.fetchone()[0]
    if max_id == None:
        max_id = 0
    
    grid_res,roa = tune_geo(elnot)
    #### insert sites
    for row in sites:
        site_incr = row[0]
        site_id = max_id + site_incr
        lat = round(row[1],3)
        lon = round(row[2],3)
        sma = round(row[3],3)
        smi = round(row[4],3)
        orient = round(row[5],1)
        heard_ct =sum(1 for i in site_ints if i[0] == site_incr)
        cur.execute("""insert into site_geos(elnot,site_id,lat,lon,sma,smi,orient,heard_count,roa) values(%s,%s,%s,%s,%s,%s,%s,%s,%s)""",
                    (elnot,site_id,lat,lon,sma,smi,orient,heard_ct,roa))
        
    #### insert site_ints
    for row in site_ints:
        site_id = max_id + row[0]
        int_id = row[1]
        cur.execute("insert into site_ints(site_id,intercept_id) values(%s,%s)",(site_id,int_id))
        
    #### Clean up residue
    int_list = rip(site_ints,1)
    if len(int_list) < max_db_list_len:
        if len(int_list)>1:
            int_str = str(tuple(int_list))
        else:
                int_str = "(" + str(int_list[0]) + ")"
        cur.execute("delete from residue_ints where r_type = 'G' and intercept_id in " + int_str)
    else:
        idx = 0
        ints_len = len(int_list)
        while idx < ints_len:
            pad = min(max_db_list_len,ints_len-idx)
            int_list2 = tuple(int_list[idx:idx+pad])
            if len(int_list2) == 1:
                int_str = "(" + str(int_list2[0]) + ")"
            else:
                int_str = str(int_list2)
            cur.execute("delete from residue_ints where r_type = 'G' and intercept_id in " + int_str)
            idx = idx + pad
        
    
    cur.execute("commit")
    cur.close()
    db.close()

    
    
def log_event(msg):
    db = db_con()
    cur=db.cursor()
    cur.execute("insert into wob_event_log(event_date, event_txt) values(now(),%s)",(msg,))
    cur.execute("commit")
    cur.close()
    db.close()


def pe_stats():
    db = db_con()
    cur=db.cursor()
    
    print()
    cur.execute("select count(1) from " + fq_idb + "intercepts")
    print("Total Intercepts: ", cur.fetchone()[0])
    
    print()
    print("***** Parametric Stats *****")
    
    cur.execute("select count(1) from parameter_groups")
    print("Parameter Groups: ", cur.fetchone()[0])
    
    
    cur.execute("select count(distinct intercept_id) from sequence_ints")
    print("Associated Intercepts (distinct): ", cur.fetchone()[0])
    
    
    cur.execute("select count(1) from sequence_ints")
    print("Associated Intercepts (including duplicates): ", cur.fetchone()[0])
    
    cur.execute("select count(1) from residue_ints where r_type='P'")
    print("Parametric Residue: ", cur.fetchone()[0])
    
    cur.execute("select count(1) from comparison_queue")
    print("Comparison_queue: ", cur.fetchone()[0])
    print()
    print("***** Geo Stats *****")
    
    cur.execute("select count(1) from site_geos")
    print("Sites: ", cur.fetchone()[0])
    
    cur.execute("select count(1) from site_ints")
    print("Site intercepts: ", cur.fetchone()[0])
    
    cur.execute("select count(1) from residue_ints where r_type='G'")
    print("Geo Residue: ", cur.fetchone()[0])
    
    print()
    print("--- End PE Stats ---")
    
    cur.close()
    db.close()
    



