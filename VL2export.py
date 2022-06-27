#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 06:44:36 2022

@author: russelmiller
"""
from VL2config import fq_refdb    # schema qualifier cnf schema
from VL2db_utils import db_con, get_tuning,config
from VL2shared_utils import trunc2
import psycopg2
#from VL2credentials import refdb_user_name,refdb_password,refdb_host,refdb_port,refdb_name


def refdb_con():
    config.read('db.config')
    
    user_name = config['refdb']['user']
    user_password = config['refdb']['password']
    db_host = config['refdb']['host']
    db_port = config['refdb']['port']
    db_name = config['refdb']['database']
    db = psycopg2.connect(user = user_name,
                          password = user_password,
                          host = db_host,
                          port = db_port,
                          database = db_name)
    return(db)



def export_pe(fq_refdb, src = 'PG', elnot_in = 0):
    src_db = db_con()
    cur = src_db.cursor()
    
    
    print()
    print("Exporting to ",fq_refdb,"....")
    
    
    if 'P' in src.upper():
        wipe_refmodes(fq_refdb,'P')
        
        #################################
        # 1. Export Global Clstrs & Scans
        #################################
        
        if elnot_in == 0:
            cur.execute("select clstr_id,elnot,mod_type,parm,parm_min,parm_max,bin_counts from global_parm_clstrs order by elnot,clstr_id")
            global_clstrs = cur.fetchall()
            cur.execute("select elnot,mod_type,scan_type from global_scans")
            global_scans = cur.fetchall()
            cur.execute("select distinct elnot,mod_type from parameter_groups")
            cases = cur.fetchall()
        else:
            cur.execute("select clstr_id,elnot,mod_type,parm,parm_min,parm_max,bin_counts from global_parm_clstrs where elnot = %s order by clstr_id",(elnot_in,))
            global_clstrs = cur.fetchall()
            cur.execute("select elnot,mod_type,scan_type from global_scans where elnot=%s",(elnot_in,))
            global_scans = cur.fetchall()
            cur.execute("select distinct '" +elnot_in+ "',mod_type from global_parm_clstrs where elnot=%s",(elnot_in,))
            cases = cur.fetchall()
            print(cases)
        
        export_global(cases,global_clstrs,global_scans,fq_refdb)
        
        
        #########################################
        # 2. Export parameter_groups (elintmodes)
        #########################################
        
        if elnot_in == 0:
            cur.execute("""select a.group_id,a.elnot,a.mod_type,b.seq_id from parameter_groups a, sequences b 
                           where a.group_id = b.group_id order by a.group_id""")
        else:
            cur.execute("""select a.group_id,a.elnot,a.mod_type,b.seq_id from parameter_groups a, sequences b 
                           where a.group_id=b.group_id and elnot = %s order by a.group_id""",(elnot_in,))
        groups = cur.fetchall()
        
        #### 2a For each group retrive sequence_elements and local clusters and export the group
        for group in groups:
            group_id = group[0]
            #elnot = group[1]
            #mod_type = group[2]
            seq_id = group[3]
            
            cur.execute("select clstr_id,position from sequence_elements where seq_id = %s",(seq_id,))
            sequence_elements = cur.fetchall()
            
            cur.execute("select clstr_id,parm,parm_min,parm_max,1 from local_parm_clstrs where group_id=%s",(group_id,))
            local_clstrs = cur.fetchall()
            
            # append global pri clstr data since it is duplicated in elintmodes
            cur.execute("""select a.clstr_id,'PRI',a.parm_min,a.parm_max,b.position from global_parm_clstrs a, sequence_elements b,
                           sequences c where a.clstr_id = b.clstr_id and b.seq_id = c.seq_id and c.group_id = %s""",(group_id,))
            local_pri = cur.fetchall()
            local_clstrs = local_clstrs + local_pri
            
            cur.execute("select scan_type from local_scans where group_id=%s",(group_id,))
            local_scans = cur.fetchall()
            
            export_group(group,sequence_elements,local_clstrs,local_scans,fq_refdb)
            
        ################################################
        # 3. Export sequence_ints (elintmode_intercepts)
        ################################################
        
        cur.execute("""select a.group_id, b.intercept_id from sequences a, sequence_ints b where a.seq_id=b.seq_id """)
        mode_ints = cur.fetchall()
        export_parm_ints(mode_ints,fq_refdb)
            
    cur.close()
    src_db.close()
    
def wipe_refmodes(fq_refdb,src = 'PG'):
    #Note: deletes ALL data from refmodes, even if only inserting a single elnot
    refdb = refdb_con()
    refcur = refdb.cursor()
    
    if 'P' in src.upper():
        refcur.execute("delete from " +fq_refdb+ "global_rf")
        refcur.execute("delete from " +fq_refdb+ "global_pri")
        refcur.execute("delete from " +fq_refdb+ "global_pd")
        refcur.execute("delete from " +fq_refdb+ "global_sp")
        refcur.execute("delete from " +fq_refdb+ "global_ir")
        
        refcur.execute("delete from " +fq_refdb+ "global_rf_bins")
        refcur.execute("delete from " +fq_refdb+ "global_pri_bins")
        refcur.execute("delete from " +fq_refdb+ "global_pd_bins")
        refcur.execute("delete from " +fq_refdb+ "global_sp_bins")
        refcur.execute("delete from " +fq_refdb+ "global_ir_bins")
        
        refcur.execute("delete from " +fq_refdb+ "globalmode")
        refcur.execute("delete from " +fq_refdb+ "globalmode_scantype")
        refcur.execute("delete from " +fq_refdb+ "elintmode")
        refcur.execute("delete from " +fq_refdb+ "elintmode_firing_order")
        refcur.execute("delete from " +fq_refdb+ "elintmode_intercepts")
        refcur.execute("delete from " +fq_refdb+ "elintmode_scantype")
        #refcur.execute("delete from " +fq_refdb+ "firingorder")
        refcur.execute("delete from " +fq_refdb+ "ir")
        refcur.execute("delete from " +fq_refdb+ "pd")
        refcur.execute("delete from " +fq_refdb+ "pri")
        refcur.execute("delete from " +fq_refdb+ "rf")
        refcur.execute("delete from " +fq_refdb+ "sp")
        #refcur.execute("delete from ")
        #refcur.execute("delete from ")
        
        
        refcur.execute("commit")
    #refcur.execute("delete from ")
    refcur.close()
    refdb.close()
    

    
def export_global(cases,global_clstrs,global_scans,fq_refdb):
    # For each case in cases (elnot, mod_type):
    #    export global clstrs (RF,PRI,PD,SP,IR): global_parm_clstrs --> global_rf, global_pri,global_pd,global_sp & global_ir
    #    populate bin structures: bin_counts --> global_rf_bins, global_pri_bins,global_pd_bins,global_sp_bins & global_ir_bins
    
    #gmodes = []
    refdb=refdb_con()
    refcur = refdb.cursor()
    
    for case in cases:
        elnot = case[0]
        mod_type = case[1]
        cwf = get_cwf(mod_type)
        if mod_type == 'B':
            iscw=True
        else:
            iscw=False
        
        # 1. Create & insert a "globalmode" for each (elnot,mod_type)
        #   - Assume there is a one to-one between (elnot,mod_type) and globalmode.  
        #   - There is no "globalmode" equivalent in VL2, so just using a one-up indx for globalmode id
        
        refcur.execute("select max(id) from " +fq_refdb+ "globalmode")
        globalmode_id = refcur.fetchone()[0]
        if globalmode_id == None:
            globalmode_id=1
        else:
            globalmode_id +=1
        
        #capture[elnot,mod_type,globalmode_id for later use
        #gmodes.append([elnot,mod_type,globalmodeid])
        
        refcur.execute("insert into " +fq_refdb+ "globalmode(id,elnot,iscw,primodtype,computedwaveform) Values(%s,%s,%s,%s,%s)",
                                              (globalmode_id,elnot,iscw,mod_type,cwf))
        
        # 2. insert global clstrs (includes bins)
        for parm in ['RF','PRI','PD','SP','IR']:
            parm_clstrs = (i for i in global_clstrs if i[1]==elnot and i[2] == mod_type and i[3] == parm)
            
            units = get_units(parm)
            incr,horizon,p = get_tuning(elnot,parm,'WOB')
            
            for clstr in parm_clstrs:
                clstr_id = clstr[0]
                parm_min = clstr[4]
                parm_max = clstr[5]
                bins = clstr[6]
                num_bins = len(bins)
                bin_min = trunc2(float(parm_min) - horizon, incr)
                bin_max = trunc2(float(parm_max) + horizon, incr)  # Note: this is the last lower bin edge, i.e., vals up to bin_max+incr fit inside the bins
                db_clstr_str = "insert into " +fq_refdb+ "global_" +parm.lower()+ "(id,typicalmin,typicalmax,units,num_bins,binmax,binmin,global_mode_id) Values(%s,%s,%s,%s,%s,%s,%s,%s)"
                refcur.execute(db_clstr_str,(clstr_id,parm_min,parm_max,units,num_bins,bin_max,bin_min,globalmode_id))
                
                # also need to update bin structures
                db_bin_str = "insert into " +fq_refdb+ "global_" +parm.lower()+ "_bins(global_" +parm.lower()+ "_id,bins_order,bins) values(%s,%s,%s)"
                bin_order = 0
                for parm_bin in bins:
                    refcur.execute(db_bin_str,(clstr_id,bin_order,parm_bin))
                    bin_order +=1
            refcur.execute("commit")

        # 3. Insert "global_mode_scan_types
        scan_types = (i[2] for i in global_scans if i[0]==elnot and i[1]==mod_type)
        for st in scan_types:
            refcur.execute("insert into " +fq_refdb+ "globalmode_scantype(globalmode_id,scan_type) values(%s,%s)",(globalmode_id,st))
            
    
    refcur.execute("commit")
    refcur.close()
    refdb.close()

    print("global export successful")
    #return gmodes

def export_group(group,sequence_elements,local_clstrs,local_scans,fq_refdb):
    
    refdb=refdb_con()
    refcur = refdb.cursor()
    
    group_id = group[0]
    elnot = group[1]
    mod_type = group[2]
    cwf = get_cwf(mod_type)
    #seq_id = group[3]
    
    # query db for globalmode_id here?
    
    # 1. insert elintmode
    if mod_type == 'B':
        iscw = True
    else:
        iscw = False
    print("exporting group",group_id)
    refcur.execute("insert into " +fq_refdb+ "elintmode(id,elnot,iscw,primodtype,reference,new_mode,computedwaveform) values(%s,%s,%s,%s,'VL2',True,%s)",
                                         (group_id,elnot,iscw,mod_type,cwf))
    
    # 2. insert sequence_elements
    for element in sequence_elements:
        clstr_id = element[0]
        position = element[1]
        print("inserting sequence_element",position)
        refcur.execute("insert into " +fq_refdb+ "elintmode_firing_order(elint_mode_id,firing_order_order,firing_order_id) values(%s,%s,%s)",(group_id,position,clstr_id))
        
    # 3. insert local clstrs
    for parm in ['RF','PRI','PD','SP','IR']:
        print("inserting local",parm,"clstrs")
        clstrs = (i for i in local_clstrs if i[1]==parm)
        db_insert_str = "insert into " +fq_refdb +parm.lower()+"(id,typicalmin,typicalmax,units,elint_mode_id,global_" +parm.lower()+ "_id," +parm.lower()+"_order)"
        db_insert_str = db_insert_str + " values(%s,%s,%s,%s,%s,%s,%s)"
        units = get_units(parm)
        refcur.execute("select max(id) from " +fq_refdb +parm.lower())
        next_id = refcur.fetchone()[0]
        if next_id == None:
            next_id = 1
        else:
            next_id +=1
        
        for clstr in clstrs:
            clstr_id = clstr[0]
            parm_min = clstr[2]
            parm_max = clstr[3]
            clstr_order = clstr[4]
            refcur.execute(db_insert_str,(next_id,parm_min,parm_max,units,group_id,clstr_id,clstr_order))
            next_id +=1
        
    
    # 4. insert local scans
    for st in local_scans:
        refcur.execute("insert into " +fq_refdb+ "elintmode_scantype(elintmode_id,scan_type) values(%s,%s)",(group_id,st))
        
    refcur.execute("commit")
    refcur.close()
    refdb.close()


def export_parm_ints(mode_ints,fq_refdb):
    refdb =refdb_con()
    refcur = refdb.cursor()
    i = 0
    for row in mode_ints:
        mode_id = row[0]
        int_id = row[1]
        refcur.execute("insert into " +fq_refdb+ "elintmode_intercepts(elint_mode_id,intercepts_intercept_id) values(%s,%s)",(mode_id,int_id))
        i +=1
        if i>100:
            refcur.execute("commit")
            i=0
    refcur.execute("commit")
    refcur.close()
    refdb.close()    
    
def get_units(parm):
    if parm.upper() == 'RF':
        units = 'MHz'
    elif parm.upper() in ['PRI','PD']:
        units = 'usec'
    elif parm.upper() == 'SP':
        units = 'sec'
    elif parm.upper() == 'IR':
        units = 'Hz'
    return units

def get_cwf(mod_type):
    # retrieve the proper text for the computed waveform field
    if mod_type.upper() =='B':
        cwf = 'CW'
    elif mod_type.upper() == 'D':
        cwf = 'CONSTANT'
    elif mod_type.upper() == 'M':
        cwf = 'STAGGER'
    elif mod_type.upper() =='S':
        cwf = 'DWELL'
    elif mod_type.upper() == 'J':
        cwf = 'JITTER'
    else:
        cwf = 'UNKNOWN'
    return cwf
            
            
def make_refdb():    
    db = refdb_con()
    cur = db.cursor()
    
    schema = fq_refdb
    if schema != '':
        schema = schema + '.'

        
    ################################### PRIMARY IPE TABLES:  #####################################
    #   Template:  cur.execute("""Create table if not exists """)
   
    cur.execute("""Create table if not exists """ +schema+ """globalmode (id int8 not null, 
                                                          computedwaveform varchar(255), 
                                                          elnot varchar(255), iscw boolean, 
                                                          mode_number varchar(255), 
                                                          primodtype varchar(255), 
                                                          priwaveform varchar(255), 
                                                          radarfunction varchar(255), 
                                                          reference varchar(255), 
                                                          rfmodtype varchar(255), 
                                                          primary key (id))""")
    
    
    cur.execute("""Create table if not exists """ +schema+ """global_rf (id int8 not null, 
                                                         extrememax float8, 
                                                         extrememin float8, 
                                                         mostprobable float8, 
                                                         name varchar(255), 
                                                         typicalmax numeric(19, 2), 
                                                         typicalmin numeric(19, 2), 
                                                         units varchar(255), 
                                                         binmax numeric(19, 2), 
                                                         binmin numeric(19, 2), 
                                                         horizon numeric(19, 2), 
                                                         increment numeric(19, 2), 
                                                         maxslope int4, 
                                                         maxslopebin int4, 
                                                         num_bins int4 not null, 
                                                         numentries int4, 
                                                         percentile numeric(19, 2), 
                                                         threshold numeric(19, 2), 
                                                         total_area numeric(19, 2), 
                                                         global_mode_id int8, 
                                                         primary key (id))""")
   
    cur.execute("""Create table if not exists """ +schema+ """global_pri (id int8 not null, 
                                                          extrememax float8, 
                                                          extrememin float8, 
                                                          mostprobable float8, 
                                                          name varchar(255), 
                                                          typicalmax numeric(19, 2), 
                                                          typicalmin numeric(19, 2), 
                                                          units varchar(255), 
                                                          binmax numeric(19, 2), 
                                                          binmin numeric(19, 2), 
                                                          horizon numeric(19, 2), 
                                                          increment numeric(19, 2), 
                                                          maxslope int4, 
                                                          maxslopebin int4, 
                                                          num_bins int4 not null, 
                                                          numentries int4, 
                                                          percentile numeric(19, 2), 
                                                          threshold numeric(19, 2), 
                                                          total_area numeric(19, 2), 
                                                          global_mode_id int8, 
                                                          pri_order int4, 
                                                          primary key (id))""")
   
    
   
    cur.execute("""Create table if not exists """ +schema+ """global_pd (id int8 not null, 
                                                         extrememax float8, 
                                                         extrememin float8, 
                                                         mostprobable float8, 
                                                         name varchar(255), 
                                                         typicalmax numeric(19, 2), 
                                                         typicalmin numeric(19, 2), 
                                                         units varchar(255), 
                                                         binmax numeric(19, 2), 
                                                         binmin numeric(19, 2), 
                                                         horizon numeric(19, 2), 
                                                         increment numeric(19, 2), 
                                                         maxslope int4, 
                                                         maxslopebin int4, 
                                                         num_bins int4 not null, 
                                                         numentries int4, 
                                                         percentile numeric(19, 2), 
                                                         threshold numeric(19, 2), 
                                                         total_area numeric(19, 2), 
                                                         global_mode_id int8, 
                                                         primary key (id))""")
    
    cur.execute("""Create table if not exists """ +schema+ """global_sp (id int8 not null, 
                                                         extrememax float8, 
                                                         extrememin float8, 
                                                         mostprobable float8, 
                                                         name varchar(255), 
                                                         typicalmax numeric(19, 2), 
                                                         typicalmin numeric(19, 2), 
                                                         units varchar(255), 
                                                         binmax numeric(19, 2), 
                                                         binmin numeric(19, 2), 
                                                         horizon numeric(19, 2), 
                                                         increment numeric(19, 2), 
                                                         maxslope int4, 
                                                         maxslopebin int4, 
                                                         num_bins int4 not null, 
                                                         numentries int4, 
                                                         percentile numeric(19, 2), 
                                                         threshold numeric(19, 2), 
                                                         total_area numeric(19, 2), 
                                                         global_mode_id int8, 
                                                         primary key (id))""")
    
    cur.execute("""Create table if not exists """ +schema+ """global_ir (id int8 not null, 
                                                         extrememax float8, 
                                                         extrememin float8, 
                                                         mostprobable float8, 
                                                         name varchar(255), 
                                                         typicalmax numeric(19, 2), 
                                                         typicalmin numeric(19, 2), 
                                                         units varchar(255), 
                                                         binmax numeric(19, 2), 
                                                         binmin numeric(19, 2), 
                                                         horizon numeric(19, 2), 
                                                         increment numeric(19, 2), 
                                                         maxslope int4, 
                                                         maxslopebin int4, 
                                                         num_bins int4 not null, 
                                                         numentries int4, 
                                                         percentile numeric(19, 2), 
                                                         threshold numeric(19, 2), 
                                                         total_area numeric(19, 2), 
                                                         global_mode_id int8, 
                                                         primary key (id))""")
    
    cur.execute("""Create table if not exists """ +schema+ """global_rf_bins (global_rf_id int8 not null, 
                                                              bins int4, 
                                                              bins_order int4 not null, 
                                                              primary key (global_rf_id, bins_order))""")
    
    cur.execute("""Create table if not exists """ +schema+ """global_pri_bins (global_pri_id int8 not null, 
                                                               bins int4, 
                                                               bins_order int4 not null, 
                                                               primary key (global_pri_id, bins_order))""")
    
    cur.execute("""Create table if not exists """ +schema+ """global_pd_bins (global_pd_id int8 not null, 
                                                              bins int4, 
                                                              bins_order int4 not null, 
                                                              primary key (global_pd_id, bins_order))""")
    
    cur.execute("""Create table if not exists """ +schema+ """global_sp_bins (global_sp_id int8 not null, 
                                                              bins int4, 
                                                              bins_order int4 not null, 
                                                              primary key (global_sp_id, bins_order))""")
    
    cur.execute("""Create table if not exists """ +schema+ """global_ir_bins (global_ir_id int8 not null, 
                                                              bins int4, 
                                                              bins_order int4 not null, 
                                                              primary key (global_ir_id, bins_order))""")
    
    cur.execute("""Create table if not exists """ +schema+ """globalmode_scantype (globalmode_id int8 not null, 
                                                                   scan_type varchar(255))""")
    
    cur.execute("""Create table if not exists """ +schema+ """elintmode (id int8 not null, 
                                                         computedwaveform varchar(255), 
                                                         elnot varchar(255), 
                                                         iscw boolean, 
                                                         mode_number varchar(255), 
                                                         primodtype varchar(255), 
                                                         priwaveform varchar(255), 
                                                         radarfunction varchar(255), 
                                                         reference varchar(255), 
                                                         rfmodtype varchar(255), 
                                                         merged boolean, 
                                                         mode_identifier varchar(255), 
                                                         new_mode boolean, 
                                                         emitter_id int8, 
                                                         mode_report_id int8, 
                                                         mode_order int4, 
                                                         source_modes_order int4, 
                                                         primary key (id))""")

    cur.execute("""Create table if not exists """ +schema+ """elintmode_firing_order (elint_mode_id int8 not null,
                                                                      firing_order_order int4 not null,
                                                                      firing_order_id int8 not null,
                                                                      primary key (elint_mode_id, firing_order_order))""")
    
    #cur.execute("""Create table if not exists firingorder (element int8 not null, 
    #                                                       firing_order_id int8 not null, 
    #                                                       firing_order_order int4 not null, 
    #                                                       primary key (element, firing_order_order))""")
    
    cur.execute("""Create table if not exists """ +schema+ """pd (id int8 not null, 
                                                  extrememax float8, 
                                                  extrememin float8, 
                                                  mostprobable float8, 
                                                  name varchar(255), 
                                                  typicalmax numeric(19, 2), 
                                                  typicalmin numeric(19, 2), 
                                                  units varchar(255), 
                                                  elint_mode_id int8, 
                                                  global_pd_id int8,
                                                  pd_order int4,
                                                  primary key (id))""")
    
    cur.execute("""Create table if not exists """ +schema+ """pri (id int8 not null, 
                                                   extrememax float8, 
                                                   extrememin float8, 
                                                   mostprobable float8, 
                                                   name varchar(255), 
                                                   typicalmax numeric(19, 2), 
                                                   typicalmin numeric(19, 2), 
                                                   units varchar(255), 
                                                   clockperiod_max float8, 
                                                   clockperiod_min float8, 
                                                   clockperiod_nominal float8, 
                                                   countdown_max int4, 
                                                   countdown_min int4, 
                                                   countdown_nomninal int4, 
                                                   is_crystalcontrolled boolean, 
                                                   mod_type varchar(255), 
                                                   elint_mode_id int8, 
                                                   global_pri_id int8, 
                                                   pri_order int4, 
                                                   primary key (id))""")
    
    cur.execute("""Create table if not exists """ +schema+ """rf (id int8 not null, 
                                                  extrememax float8, 
                                                  extrememin float8, 
                                                  mostprobable float8, 
                                                  name varchar(255), 
                                                  typicalmax numeric(19, 2), 
                                                  typicalmin numeric(19, 2), 
                                                  units varchar(255), 
                                                  mod_type varchar(255), 
                                                  elint_mode_id int8, 
                                                  global_rf_id int8, 
                                                  rf_order int4,
                                                  primary key (id))""")
    
    cur.execute("""Create table if not exists """ +schema+ """sp (id int8 not null, 
                                                  extrememax float8, 
                                                  extrememin float8, 
                                                  mostprobable float8, 
                                                  name varchar(255), 
                                                  typicalmax numeric(19, 2), 
                                                  typicalmin numeric(19, 2), 
                                                  units varchar(255), 
                                                  elint_mode_id int8, 
                                                  global_sp_id int8,
                                                  sp_order int4,
                                                  primary key (id))""")
    
    cur.execute("""Create table if not exists """ +schema+ """ir (id int8 not null, 
                                                  extrememax float8, 
                                                  extrememin float8, 
                                                  mostprobable float8, 
                                                  name varchar(255), 
                                                  typicalmax numeric(19, 2), 
                                                  typicalmin numeric(19, 2), 
                                                  units varchar(255), 
                                                  elint_mode_id int8, 
                                                  global_ir_id int8,
                                                  ir_order int4,
                                                  primary key (id))""")
    
    cur.execute("""Create table if not exists """ +schema+ """elintmode_intercepts (elint_mode_id int8 not null, 
                                                                    intercepts_intercept_id int8 not null)""")
    
    cur.execute("""Create table if not exists """ +schema+ """elintmode_scantype (elintmode_id int8 not null, 
                                                                  scan_type varchar(255))""")
    
    #cur.execute("""Create table if not exists """ +schema+ """ """)
    cur.execute("Commit")
    cur.close()
    db.close()
    
            
            
            
            
            
