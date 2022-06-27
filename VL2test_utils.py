#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 16:20:31 2021
@author: russelmiller


int_bldr(num_ints, delete_prev = True, insert_cq = True, seeds = [])
    - create randomized intercept data for testing

int_clstr_histo(elnot,parm,mod_type,clstr_id)
    - plot histogram of data in the vicinity of a specified cluster

int_histo(elnot,parm,mod_type,parm_min,parm_max)
    - plot histogram based on all available intercepts
    
int_mode_histo(mode_id,parm,parm_min = None, parm_max = None)
    - plot histograms for intercepts attached to a specified mode

res_histo(elnot,parm,mod_type,parm_min,parm_max):
    - plot histograms based on date in residue

"""

from VL2db_utils import db_con, get_tuning,delete_env
import numpy as np
import matplotlib.pyplot as plt
from VL2shared_utils import trunc2, rip
from VL2config import fq_idb    # ONLY used for cluster histograns, and NOT int_bldr (int generation in local schema only)

def int_bldr(num_ints, delete_prev = True,insert_cq = True, seeds_in = []):
    db = db_con()
    cur = db.cursor()
    from numpy.random import normal as N
    import random
    
    system_random=random.SystemRandom()
    
    max_elems = 2       # maximum number of positions M/S/J
    min_elems = 2       # default = 2...set min & max equal to generate sequences all the same size for M/S/J
    geo_dev = 0.05      # Max magnitude for random geo variation (decimal degrees)

    if seeds_in != []:
        seeds = seeds_in
    else:
        seeds = [
            ['ELNOT',['BBBBB']],
            ['MT',['M','D','S']],
            ['ST',['A','Z']],
            ['RF',[3000,4000],[1,200]],
            ['RF',[5000],[10,500]],
            ['PRI',[100],[1.5,200]],
            ['PRI',[200],[2,200]],
            #['PD',[1.5],[0.15,100]],
            #['PD',[5,10],[0.1,100]],
            ['SP',[5,10,15],[0.1,100]],
            ['IR',[20],[1,200]],
            #['GEO', [[45.0,120.0],[47.0,125.0],[50.0,130.0],[50.0,120.0],[46.0,140.0]],[.05]]
            #['GEO', [[29.0,-98.0],[39.0,-105.0]],[.05]]
            ['GEO', [[29.0,-98.0],[39.0,-104.0]], [0.05]]
            ]
    #print(seeds)
    elnot_seeds = []
    mt_seeds = []
    st_seeds = []
    geo_seeds=[]
    rf_dat = []
    pri_dat = []
    pd_dat = []
    sp_dat = []
    ir_dat = []
    for seed in seeds:
        parm = seed[0].upper()
        if parm in ['ELNOT','MT','ST','GEO']:
            new_vals = seed[1]
            if parm == 'ELNOT':
                elnot_seeds = elnot_seeds + new_vals
            elif parm == 'MT':
                mt_seeds = mt_seeds + new_vals
            elif parm == 'ST':
                st_seeds = st_seeds + new_vals
            elif parm =='GEO':
                geo_seeds = seed[1]
                geo_dev = seed[2][0]
                
        elif parm != 'GEO':
            std = seed[2][0]
            cnt = seed[2][1]
            for value in seed[1]:
                dat = N(value,std, size=cnt)
                if parm =='RF':
                    rf_dat = np.append(rf_dat,dat)
                elif parm == 'PRI':
                    pri_dat = np.append(pri_dat,dat)
                elif parm == 'PD':
                    pd_dat = np.append(pd_dat,dat)
                elif parm == 'SP':
                    sp_dat = np.append(sp_dat,dat)
                elif parm == 'IR':
                    ir_dat = np.append(ir_dat,dat)
        
            


    if delete_prev == True:
        delete_env(0)
        cur.execute("delete from intercept_pris")
        cur.execute("delete from intercepts")
        
        next_id=1024
    else:
        cur.execute("select max(intercept_id) from intercepts")
        next_id = cur.fetchone()[0]
        if next_id == None:
            next_id =1024
        else:
            next_id +=1
    
    
    for ints in range(next_id,next_id+num_ints):
        int_id = ints

        # Generate/Insert Intercepts
        elnot = system_random.choice(elnot_seeds)
        mt_val = system_random.choice(mt_seeds)
        
        rf_val = round(system_random.choice(rf_dat),3)
        if rf_val < 0:
            continue
        if len(pd_dat) > 0:
            pd_val = round(system_random.choice(pd_dat),3)
        else:
            pd_val = -1
        if len(sp_dat) > 0:
            sp_val = round(system_random.choice(sp_dat),3)
        else:
            sp_val = -1
        if len(ir_dat) > 0:
            ir_val = round(system_random.choice(ir_dat),3)
        else:
            ir_val = -1
        if len(st_seeds) > 0:
            st_val = system_random.choice(st_seeds)
        else:
            st_val = '-'
        
        
        if len(geo_seeds)>0:
            sma = system_random.choice(range(200,500))/100
            smi = system_random.choice(range(10,200))/100
            orient = system_random.choice(range(0,1800))/10
            coords = system_random.choice(geo_seeds)
            lat = coords[0] +geo_dev * system_random.choice(range(-100,100))/100
            lon = coords[1] +geo_dev * system_random.choice(range(-100,100))/100
        else:
            sma=-1
            smi=-1
            orient=-1
            lat=-1
            lon=-1
        
      
        cur.execute("""insert into intercepts(intercept_id,elnot,mod_type,rf,pd,sp,ir,scan_type,latitude,longitude,orientation,area_dist_maj,
                    area_dist_min,wrangler_id, is_emitter,time_process) 
                    values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, 'Y',current_date)""",(int_id,elnot,mt_val,rf_val,pd_val,sp_val,ir_val,st_val,lat,lon,orient,sma,smi,int_id))
        
        # Generate/Insert PRIs:
        
        # Determine sequence length
        if mt_val =='D':
            slen = 1
        else:
            slen = system_random.choice(range(min_elems,max_elems+1))
            
        # Generate/Insert PRIs
        for i in range(0,slen):
            position = i + 1
            
            pri_val = round(system_random.choice(pri_dat),3)
            cur.execute("insert into intercept_pris(intercept_id,pri_value,pri_number) values(%s,%s,%s)",(int_id,pri_val,position))
            
    if insert_cq == True:
        cur.execute("insert into comparison_queue(intercept_id) select intercept_id from intercepts where intercept_id >= %s",(next_id,))

        
    cur.execute("commit")
    cur.close()
    db.close()

def int_clstr_histo(elnot,parm,mod_type,clstr_id):
    # go 2 horizons either side of specified clstr
    db = db_con()
    cur = db.cursor()
    
    incr,horizon,p = get_tuning(elnot,parm)
    
    if parm == 'RF':
        cur.execute("select rf_min,rf_max from wob_rf_clstrs where elnot=%s and mod_type=%s and clstr_id=%s",
                                (elnot,mod_type,clstr_id))
        dat = cur.fetchone()
        if dat == None:
            print("specified cluster does not exist")
            return
        clstr_min = float(dat[0])
        clstr_max = float(dat[1])
        cur.execute("select rf from " + fq_idb + "intercepts where elnot=%s and adj_mt(mod_type)=%s and rf between %s and %s",
                    (elnot,mod_type, clstr_min-horizon,clstr_max+horizon))
    elif parm == 'PRI':
        cur.execute("select pri_min,pri_max from wob_pri_clstrs where elnot=%s and mod_type=%s and clstr_id=%s",
                                (elnot,mod_type,clstr_id))
        dat = cur.fetchone()
        if dat == None:
            print("specified cluster does not exist")
            return
        clstr_min = float(dat[0])
        clstr_max = float(dat[1])
        cur.execute("""select pri_value from """ + fq_idb + """intercepts a, """ + fq_idb + """intercept_pris b where a.intercept_id = b.intercept_id and elnot=%s 
                    and adj_mt(a.mod_type)=%s and pri_value between %s and %s""",
                    (elnot,mod_type, clstr_min-horizon,clstr_max+horizon))
    elif parm == 'PD':
        cur.execute("select parm_min,parm_max from wob_mode_clstrs a, wob_modes b where a.mode_id = b.mode_id and b.elnot=%s and b.mod_type =%s and a.parm = %s and a.clstr_id=%s",
                                (elnot,mod_type,parm,clstr_id))
        dat = cur.fetchone()
        if dat == None:
            print("specified cluster does not exist")
            return
        clstr_min = float(dat[0])
        clstr_max = float(dat[1])
        cur.execute("select pd from " + fq_idb + "intercepts where elnot=%s and adj_mt(mod_type)=%s and pd between %s and %s",
                    (elnot,mod_type, clstr_min-horizon,clstr_max+horizon))
    elif parm == 'SP':
        cur.execute("select parm_min,parm_max from wob_mode_clstrs a, wob_modes b where a.mode_id = b.mode_id and b.elnot=%s and b.mod_type =%s and a.parm = %s and a.clstr_id=%s",
                                (elnot,mod_type,parm,clstr_id))
        dat = cur.fetchone()
        if dat == None:
            print("specified cluster does not exist")
            return
        clstr_min = float(dat[0])
        clstr_max = float(dat[1])
        cur.execute("select sp from " + fq_idb + "intercepts where elnot=%s and adj_mt(mod_type)=%s and sp between %s and %s",
                    (elnot,mod_type, clstr_min-horizon,clstr_max+horizon))
    elif parm == 'IR':
        cur.execute("select parm_min,parm_max from wob_mode_clstrs a, wob_modes b where a.mode_id = b.mode_id and b.elnot=%s and b.mod_type =%s and a.parm = %s and a.clstr_id=%s",
                                (elnot,mod_type,parm,clstr_id))
        
        dat = cur.fetchone()
        if dat == None:
            print("specified cluster does not exist")
            return
        clstr_min = float(dat[0])
        clstr_max = float(dat[1])
        cur.execute("select ir from " + fq_idb + "intercepts where elnot=%s and adj_mt(mod_type)=%s and ir between %s and %s",
                    (elnot,mod_type, clstr_min-horizon,clstr_max+horizon))
    else:
        print('parm not supported: ',parm)
        return
    
    dat = cur.fetchall()
    if dat == []:
        print("No "+parm+" intercept data exists in the range "+str(clstr_min-horizon)+" - "+str(clstr_max+horizon))
        return
    else:
        print(len(dat),"values retrieved")
    
    dat = np.array(dat).astype(float)

    lo_val = clstr_min-horizon
    hi_val = clstr_max+horizon
    print("horizon,bin_min,bin_max,clstr_min,clstr_max=",horizon,lo_val,hi_val,clstr_min,clstr_max)
    n_bins = int(10 * (trunc2(hi_val,.1) - trunc2(lo_val,.1)))  # default is 0.1 bin width
    if n_bins > 1000:
        n_bins = int(n_bins/10)
    print("nbins",n_bins)
    hist, bin_edges = np.histogram(dat, bins = n_bins)

     
    plt.hist(dat, bins = bin_edges)
    x_lab = elnot + " " + parm + " MT=" + mod_type + " (p=" +str(p) + "; horizon=" + str(horizon) + ")"
    title = "clstr_id=" + str(clstr_id) + ": [" + str(clstr_min) + ", " + str(clstr_max) + ")"
    plt.title(title)
    plt.xlabel(x_lab)
    plt.show()
    
    
def int_histo(elnot,parm,mod_type,parm_min,parm_max):
    # go 2 horizons either side of specified clstr
    db = db_con()
    cur = db.cursor()
    
    incr,horizon,p = get_tuning(elnot,parm)
    
    if parm == 'RF':
        cur.execute("select rf from " + fq_idb + "intercepts where elnot=%s and adj_mt(mod_type)=%s and rf between %s and %s",
                    (elnot,mod_type, parm_min,parm_max))
    elif parm == 'PRI':
        cur.execute("""select pri_value from """ + fq_idb + """intercepts a, """ + fq_idb + """intercept_pris b where a.intercept_id = b.intercept_id and elnot=%s 
                    and adj_mt(a.mod_type)=%s and pri_value between %s and %s""",
                    (elnot,mod_type, parm_min,parm_max))
    elif parm == 'PD':
        cur.execute("select pd from " + fq_idb + "intercepts where elnot=%s and adj_mt(mod_type)=%s and pd between %s and %s",
                    (elnot,mod_type, parm_min,parm_max))
    elif parm == 'SP':
        cur.execute("select sp from " + fq_idb + "intercepts where elnot=%s and adj_mt(mod_type)=%s and sp between %s and %s",
                    (elnot,mod_type, parm_min,parm_max))
    elif parm == 'IR':
        cur.execute("select ir from " + fq_idb + "intercepts where elnot=%s and adj_mt(mod_type)=%s and ir between %s and %s",
                    (elnot,mod_type, parm_min,parm_max))
        
    dat = cur.fetchall()
    if dat == []:
        print("No "+parm+" intercept data exists in the range "+str(parm_min)+" - "+str(parm_max))
        return
    else:
        print(len(dat),"values retrieved")
    
    dat = np.array(dat).astype(float)
    lo_val = parm_min
    hi_val = parm_max
    n_bins = int(10 * (trunc2(hi_val,.1) - trunc2(lo_val,.1)))  # default is 0.1 bin width
    if n_bins > 1000:
        n_bins = int(n_bins/10)
    print("nbins",n_bins)
    hist, bin_edges = np.histogram(dat, bins = n_bins)
     
    plt.hist(dat, bins = bin_edges)
    x_lab = elnot + " " + parm + " MT=" + mod_type + ")"
    #title = "clstr_id=" + str(clstr_id) + ": [" + str(clstr_min) + ", " + str(clstr_max) + ")"
    #plt.title(title)
    plt.xlabel(x_lab)
    plt.show()



def int_mode_histo(mode_id,parm,parm_min = None, parm_max = None):

    db = db_con()
    cur = db.cursor()
    
    
    if parm == 'RF' or parm == 'RFL':
        cur.execute("""select rf from """ + fq_idb + """intercepts a, sequence_ints b, sequences c where 
                    a.intercept_id = b.intercept_id and b.seq_id = c.seq_id and c.group_id=%s""",(mode_id,))
    elif parm == 'PRI':
        cur.execute("""select pri_value from """ + fq_idb + """intercepts a, """ + fq_idb + """intercept_pris b, sequence_ints c, sequences d
                     where a.intercept_id = b.intercept_id and a.intercept_id= c.intercept_id and c.seq_id = d.seq_id 
                     and d.group_id = %s""",(mode_id,))
    elif parm == 'PD':
        cur.execute("""select pd from """ + fq_idb + """intercepts a, sequence_ints b, sequences c where 
                    a.intercept_id = b.intercept_id and b.seq_id = c.seq_id and c.group_id=%s and pd > 0 """,(mode_id,))
    elif parm == 'SP':
        cur.execute("""select sp from """ + fq_idb + """intercepts a, sequence_ints b, sequences c where 
                    a.intercept_id = b.intercept_id and b.seq_id = c.seq_id and c.group_id=%s and sp > 0""",(mode_id,))
    elif parm == 'IR':
        cur.execute("""select ir from """ + fq_idb + """intercepts a, sequence_ints b, sequences c where 
                    a.intercept_id = b.intercept_id and b.seq_id = c.seq_id and c.group_id=%s and ir > 0""",(mode_id,))
        
    dat = cur.fetchall()
    dat = rip(dat,0)
    if len(dat) == 0:
        print("No data found...exiting")
        return
    if parm_min == None:
        parm_min = min(dat)
    if parm_max == None:
        parm_max = max(dat)
    
    
    if dat == []:
        print("No data found")
        #print("No "+parm+" intercept data exists in the range "+str(parm_min)+" - "+str(parm_max))
        return
    else:
        print(len(dat),"values retrieved")
    
    dat = np.array(dat).astype(float)
    lo_val = parm_min
    hi_val = parm_max
    #print("horizon,bin_min,bin_max,clstr_min,clstr_max=",horizon,lo_val,hi_val,clstr_min,clstr_max)
    n_bins = int(10 * (trunc2(hi_val,.1) - trunc2(lo_val,.1)))  # default is 0.1 bin width
    if n_bins > 1000:
        n_bins = int(n_bins/10)
    print("nbins",n_bins)
    hist, bin_edges = np.histogram(dat, bins = n_bins)
     
    plt.hist(dat, bins = bin_edges)
    #x_lab = elnot + " " + parm + " MT=" + mod_type + ")"
    title = parm + "    mode_id=" + str(mode_id)
    plt.title(title)
    #plt.xlabel(x_lab)
    plt.show()







def res_histo(elnot,parm,mod_type,parm_min,parm_max):
    # go 2 horizons either side of specified clstr
    db = db_con()
    cur = db.cursor()
    
    incr,horizon,p = get_tuning(elnot,parm)
    
    if parm == 'RF':
        cur.execute("select rf from " + fq_idb + "intercepts a, residue_ints b where a.intercept_id = b.intercept_id and elnot=%s and adj_mt(mod_type)=%s and rf between %s and %s and r_type = 'P'",
                    (elnot,mod_type, parm_min,parm_max))
    elif parm == 'PRI':
        cur.execute("""select pri_value from """ + fq_idb + """intercepts a, """ + fq_idb + """intercept_pris b , residue_ints c where a.intercept_id = b.intercept_id and 
                    a.intercept_id = c.intercept_id and elnot=%s and adj_mt(a.mod_type)=%s and pri_value between %s and %s and r_type = 'P'""",
                    (elnot,mod_type, parm_min,parm_max))
    elif parm == 'PD':
        cur.execute("select pd from " + fq_idb + """intercepts a, residue_ints b where a,intercept_id = b.intercept_id and elnot=%s and adj_mt(mod_type)=%s 
                    and pd between %s and %s and r_type = 'P' """,
                    (elnot,mod_type, parm_min,parm_max))
    elif parm == 'SP':
        cur.execute("select sp from " + fq_idb + """intercepts a, residue_ints b where a.intercept_id = b.intercept_id and elnot=%s and adj_mt(mod_type)=%s 
                    and sp between %s and %s and r_type='P' """,
                    (elnot,mod_type, parm_min,parm_max))
    elif parm == 'IR':
        cur.execute("select ir from " + fq_idb + """intercepts a, residue_ints b where a.intercept_id = b.intercept_id and elnot=%s and adj_mt(mod_type)=%s 
                    and ir between %s and %s and r_type = 'P' """,
                    (elnot,mod_type, parm_min,parm_max))
        
    dat = cur.fetchall()
    if dat == []:
        print("No "+parm+" intercept data exists in the range "+str(parm_min)+" - "+str(parm_max))
        return
    else:
        print(len(dat),"values retrieved")
    
    dat = np.array(dat).astype(float)
    lo_val = parm_min
    hi_val = parm_max
    #print("horizon,bin_min,bin_max,clstr_min,clstr_max=",horizon,lo_val,hi_val,clstr_min,clstr_max)
    n_bins = int(10 * (trunc2(hi_val,.1) - trunc2(lo_val,.1)))  # default is 0.1 bin width
    if n_bins > 1000:
        n_bins = int(n_bins/10)
    print("nbins",n_bins)
    hist, bin_edges = np.histogram(dat, bins = n_bins)
    #return hist,bin_edges

     
    plt.hist(dat, bins = bin_edges)
    x_lab = elnot + " " + parm + " MT=" + mod_type + ")"
    #title = "clstr_id=" + str(clstr_id) + ": [" + str(clstr_min) + ", " + str(clstr_max) + ")"
    #plt.title(title)
    plt.xlabel(x_lab)
    plt.show()
