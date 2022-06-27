#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 09:58:40 2021

@author: russelmiller
"""
from VL2db_utils import db_con, get_tuning, delete_env, add_residue,log_event
from VL2config import min_baseline_ints, min_clstr_ints, fq_idb, baseline_geo, gb_iters,max_SMA

from VL2dm_utils import wobulate_local_seed, purge_conformal,cluster_data
from VL2shared_utils import rip, simple_merge
from VL2bin_utils import append_clstr_bins
from VL2geominer import proc_geo
from VL2geo_associator import associate_geo

from VL2residue_miner import db_global_parm_insert,new_residue_miner
from datetime import datetime


# Need to reconcile and determine which version to use
#from residue_miner import cluster_data


## Baseline Builder
##
## 1. iteratively cluster and stabilize RF & PRI
## 2. create & populate RF & PRI bin structures (scratch_pad.initial_bin_build)
## 3. stage data in residue
## 4. run residue_miner
##    4a. update residue miner to create & populate bins for new clstrs (scratch_pad.initial_bin_build)

def build_baseline(elnot,include_geo=True):
    start = datetime.now()
    db = db_con()
    cur = db.cursor()
    
    #============================================================
    # 1.  Delete any previous characterization and staged data
    #============================================================
    
    delete_env(elnot)
    cur.execute("""delete from comparison_queue where intercept_id in (select a.intercept_id from comparison_queue a,
                """ + fq_idb + """intercepts b where a.intercept_id=b.intercept_id and b.elnot = %s)""",(elnot,))
    cur.execute("""delete from residue_ints where intercept_id in (select intercept_id from """ + fq_idb + """intercepts 
                where elnot = %s)""",(elnot,))
    cur.execute("commit")
    
    
    
    #======================================================================
    # 2. Build global clusters for each mod_type having sufficient data
    #======================================================================
    
    cur.execute("select adj_mt(mod_type) mt, count(*) from " + fq_idb + "intercepts where elnot = %s group by mt having count(*) >= %s order by mt",
                (elnot,min_baseline_ints))
    cases = cur.fetchall()
    
    for case in cases:
        mod_type = case[0]
        
        #### Build & store mature clusters and stage individual parameter residue
        for parm in ['RF','PRI','PD','SP','IR']:
            baseline_parm(elnot,mod_type,parm)  # builds & stores mature global clstrs, and populates parm_residue
    print("Global parm clustering complete...")
    

    #=============================================
    # 3. Call residue_miner to build out modes
    #=============================================
        
    cur.execute("""insert into residue_ints(intercept_id,r_type) select intercept_id,'P' from """ + fq_idb + """intercepts where elnot = %s""",(elnot,))
    cur.execute("commit")
    new_residue_miner(elnot,False)
        

    #===========================================================
    # 4. Build Geo baseline if baseline_geo = True in VLconfig
    #===========================================================
    
    if baseline_geo and include_geo and elnot[0] not in ['A','C','K','M','O']:
        cur.execute("insert into residue_ints select intercept_id,'G' from " + fq_idb + """intercepts where elnot = %s
                       and latitude is not null and longitude is not null and area_dist_maj is not null and area_dist_min is not null
                       and area_dist_maj > 0 and area_dist_min > 0 and orientation >= 0 
                       and orientation is not null and is_emitter = 'Y' and area_dist_maj <= %s""",(elnot,max_SMA))
        cur.execute("commit")
        for i in range(0,gb_iters):
            sites,site_ints = proc_geo(elnot)
            if sites == []:
                print("Geo stabilized after",i,"iterations")
                break
            else:
                print(len(sites),"new geolocations retrieved...continuing")
                cur.execute("delete from comparison_geo_snapshot")
                cur.execute("insert into comparison_geo_snapshot select b.intercept_id from " + fq_idb + """
                            intercepts a, residue_ints b where
                            a.intercept_id = b.intercept_id and a.elnot = %s and b.r_type = 'G'""",(elnot,))
                cur.execute("commit")
                cur.execute("delete from residue_ints where r_type = 'G' and intercept_id in (select intercept_id from comparison_geo_snapshot)")
                cur.execute("commit")
                cur.execute("select count(1) from residue_ints where r_type='G'")
                rcnt =cur.fetchone()[0]
                print(rcnt,"ints in residue")
                cur.execute("select count(1) from comparison_geo_snapshot")
                rcnt =cur.fetchone()[0]
                print(rcnt,"ints in comp_geo_snapshot....calling associator")
                associate_geo()
                cur.execute("select count(1) from residue_ints where r_type='G'")
                rcnt =cur.fetchone()[0]
                print(rcnt,"ints in residue")
    
    cur.close()
    db.close()
    end = datetime.now()
    print("build_baseline start: ",start)
    print("build_baseline end: ",end)
    print("build_baseline total:",end-start)
    ### end of build_baseline()
    
    

def baseline_parm(elnot,mod_type,parm):
    db = db_con()
    cur = db.cursor()
    
    ## Retrieve data
    if parm == 'RF':
        cur.execute(""" select rf from """ + fq_idb + """intercepts where elnot = %s and adj_mt(mod_type) = %s order by rf""",(elnot,mod_type))
    elif parm == 'PD':
        cur.execute(""" select pd from """ + fq_idb + """intercepts where elnot = %s and adj_mt(mod_type) = %s and pd>0 order by pd""",(elnot,mod_type))
    elif parm == 'SP':
        cur.execute(""" select sp from """ + fq_idb + """intercepts where elnot = %s and adj_mt(mod_type) = %s and sp>0 order by sp""",(elnot,mod_type))
    elif parm == 'IR':
        cur.execute(""" select ir from """ + fq_idb + """intercepts where elnot = %s and adj_mt(mod_type) = %s and ir>0 order by ir""",(elnot,mod_type))
    elif parm == 'PRI':
        cur.execute(""" select pri_value from """ + fq_idb + """intercept_pris a, """ + fq_idb + """intercepts b where a.intercept_id = b.intercept_id and
                        b.elnot = %s and adj_mt(b.mod_type) = %s order by pri_value""",(elnot,mod_type))
    dat = rip(cur.fetchall(),0)
    
    if dat == []:
        return dat
    
    incr,horizon,thresh,p = get_tuning(elnot,parm,'BOTH')
    

    ########
    #### Extract "mature" histogram-based clusters 
    ####   - repeat clustering until no new clusters are found or there is not enough data left to process
    ####   - merge new clstrs with existing to maintain minimum separation (horizon)
    ########
    
    clstrs = cluster_data(dat,incr,horizon,thresh)
    clstrs = simple_merge(clstrs,horizon)
    residue = purge_conformal(dat,clstrs,horizon)
            
    if len(residue) >= min_clstr_ints:
        maturity = False
    else:
        maturity = True
    ebrake = 0
    while not maturity:
        ebrake +=1
        if ebrake>20:
            event_str = "####WARNING###: Failed to mature " +parm+ " clstr: " +elnot+ " " +mod_type+ " (VL2bb.baseline_parm-1)"
            print(event_str)
            log_event(event_str)
            break
        new_clstrs = cluster_data(residue,incr,horizon,thresh)
        if len(new_clstrs) == 0:
            maturity = True
        else:
            for row in new_clstrs:
                clstrs.append(row)
            clstrs = simple_merge(clstrs, horizon)
            residue = purge_conformal(dat,clstrs,horizon)
        if len(residue) < min_clstr_ints:
            maturity = True            
            
            
    #######
    #### Wobulate & stabilize new clusters (assymetrical percentile-based (ASP) algorithm)
    ########
    
    merge_flag = True
    ebrake = 0
    while merge_flag == True:
        ebrake +=1
        if ebrake>20:
            event_str = "####WARNING###: Failed to stabilize " +parm+ " clstr: " +elnot+ " " +mod_type+ " (VL2bb.baseline_parm-2)"
            print(event_str)
            log_event(event_str)
            break
        wob_clstrs = []
        for clstr in clstrs:
            clstr_min_seed = clstr[0]
            clstr_max_seed = clstr[1]
            c_min,c_max = wobulate_local_seed(dat,clstr_min_seed, clstr_max_seed,incr,horizon,p)
            # Note: wobulate_local_seed (vs wobulate_global_seed) since all data is available and no db clstrs yet exist
            wob_clstrs.append([c_min,c_max])
        clstrs = simple_merge(wob_clstrs,horizon)  # ensure clstrs are adequately separated
        if len(clstrs) == len(wob_clstrs):
            merge_flag = False

    ########
    #### Append bins and store new clstrs
    ########

    parm_inserts = append_clstr_bins(dat,clstrs,horizon,incr,p)  # parm_inserts is of form [clstr_min,clstr_max,[bins]]
    db_global_parm_insert(elnot,mod_type,parm,parm_inserts)
    
    parm_residue = purge_conformal(dat,clstrs,horizon)
    add_residue(elnot,mod_type,parm,parm_residue,incr)
    
    
    cur.execute("commit")
    cur.close()
    db.close()
    #return clstrs












