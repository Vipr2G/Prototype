#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:02:57 2021

@author: russelmiller
"""
from VL2db_utils import db_con
from VL2shared_utils import add_unique

def pretty_print(elnot = None):
    """
    Prints easily readable summary of modes in the database
    """
    db = db_con()
    cur = db.cursor()
    if elnot == None:
        cur.execute("select group_id,elnot, mod_type from parameter_groups order by group_id")
    else:
        cur.execute("select group_id,elnot, mod_type from parameter_groups where elnot = %s order by group_id",(elnot,))
    modes = cur.fetchall()
    for row in modes:
        #unique_elems = []
        mode_id = row[0]
        elnot = row[1]
        mod_type = row[2]
        
        cur2 = db.cursor()
        cur2.execute(" select count(a.intercept_id) from sequence_ints a, sequences b where a.seq_id = b.seq_id and b.group_id = %s",(mode_id,))
        mode_hc = cur2.fetchone()[0]
        cur2.close()
        
        print()
        print('-----------------------------------------------------------------------------------')
        print()

        
        p_str = "group_id: " + str(mode_id) + "     ELNOT: " + elnot + "     Mod Type: " + mt_str(mod_type) + "     HEARD COUNT: " + str(mode_hc)
        print(p_str)
        print()

        
        # get parent sequence
        cur.execute(""" select a.clstr_id, b.parm_min,b.parm_max from sequence_elements a, global_parm_clstrs b, sequences c 
                    where a.clstr_id = b.clstr_id  and a.seq_id =c.seq_id and c.group_id = %s and c.is_parent = 'T' order by a.position""",(mode_id,))
        pri_elems = cur.fetchall()
        seq,seq_clstrs = firing_order(pri_elems)
        seq_str = ''
        for elem in seq:
            if seq_str != '':
                seq_str = seq_str + '-'
            seq_str = seq_str + str(elem) 
        print("Firing Order:",seq_str)
        print()
        #print(seq,seq_clstrs)
        for i in range(0,len(seq_clstrs)):
            if i == 0:
                p_str = "PRI CLUSTERS:     " + str(seq_clstrs[i][0]) + " - " + str(seq_clstrs[i][1])
            else:
                p_str = "                  " + str(seq_clstrs[i][0]) + " - " + str(seq_clstrs[i][1])
            print(p_str)
        print()
        
        #cur.execute("""select a.parm_min, a.parm_max from global_parm_clstrs a, local_parm_clstrs b where a.clstr_id = b.clstr_id 
        #            and a.parm = 'RF' and b. group_id = %s order by a.parm_min""",(mode_id,))
        #rf_clstrs = cur.fetchall()
        
        ## Local Clusters
        
        cur.execute("select parm_min, parm_max from local_parm_clstrs where parm = 'RF' and group_id = %s order by parm_min",(mode_id,))
        rf_clstrs = cur.fetchall()
        cur.execute("select parm_min, parm_max from local_parm_clstrs where parm = 'PD' and group_id = %s order by parm_min",(mode_id,))
        pd_clstrs = cur.fetchall()
        cur.execute("select parm_min, parm_max from local_parm_clstrs where parm = 'SP' and group_id = %s order by parm_min",(mode_id,))
        sp_clstrs = cur.fetchall()
        cur.execute("select parm_min, parm_max from local_parm_clstrs where parm = 'IR' and group_id = %s order by parm_min",(mode_id,))
        ir_clstrs = cur.fetchall()
        cur.execute("select scan_type from local_scans where group_id = %s order by scan_type",(mode_id,))
        s_types = cur.fetchall()
        
        # Global Clusters
        cur.execute("""select a.parm_min,a.parm_max from global_parm_clstrs a,local_parm_clstrs b where a.clstr_id=b.clstr_id
                    and b.parm='RF' and b.group_id = %s order by a.parm_min""",(mode_id,))
        global_rf = cur.fetchall()
        
        cur.execute("""select a.parm_min,a.parm_max from global_parm_clstrs a,local_parm_clstrs b where a.clstr_id=b.clstr_id
                    and b.parm='PD' and b.group_id = %s order by a.parm_min""",(mode_id,))
        global_pd = cur.fetchall()
        
        cur.execute("""select a.parm_min,a.parm_max from global_parm_clstrs a,local_parm_clstrs b where a.clstr_id=b.clstr_id
                    and b.parm='SP' and b.group_id = %s order by a.parm_min""",(mode_id,))
        global_sp = cur.fetchall()
        
        cur.execute("""select a.parm_min,a.parm_max from global_parm_clstrs a,local_parm_clstrs b where a.clstr_id=b.clstr_id
                    and b.parm='IR' and b.group_id = %s order by a.parm_min""",(mode_id,))
        global_ir = cur.fetchall()
        l1 = len(rf_clstrs)
        #l11 = len(rfl_clstrs)
        l2 = len(pd_clstrs)
        l3 = len(sp_clstrs)
        l4 = len(ir_clstrs)
        l5 = len(s_types)
        max_len = max(l1,l2,l3,l4,l5)
        print("Clusters:  local (global)")
        for line in range(1,max_len+1):
            if line == 1:
                rf_str = 'RFs: '
                #rfl_str = '    RFs (local): '
                pd_str = '     PDs: '
                sp_str = '     SPs: '
                st_str = '     STs: '
                ir_str = '     IRs: '
            else:
                rf_str = '              '
                #rfl_str = '                 '
                pd_str = '              '
                sp_str = '          '
                st_str = '            '
                ir_str = '          '
                
            if l1 >= line:
                rf_str = rf_str + str(rf_clstrs[line-1][0]) + ' - ' + str(rf_clstrs[line-1][1])
                rf_str = rf_str + ' (' + str(global_rf[line-1][0]) + ' - ' + str(global_rf[line-1][1]) + ')'
            else:
                rf_str = '                                  '
                
            
            #if l11 >= line:
            #    rfl_str = rfl_str + str(rfl_clstrs[line-1][0]) + ' - ' + str(rfl_clstrs[line-1][1])
            #else:
            #    rfl_str = '                                '
            
            
            
            if l2 >= line:
                pd_str = pd_str + str(pd_clstrs[line-1][0]) + ' - ' + str(pd_clstrs[line-1][1])
                pd_str = pd_str + ' (' + str(global_pd[line-1][0]) + ' - ' + str(global_pd[line-1][1]) + ')'
            else:
                pd_str = '                                     '
            
            if l3 >= line:
                sp_str = sp_str + str(sp_clstrs[line-1][0]) + ' - ' + str(sp_clstrs[line-1][1])
                sp_str = sp_str + ' (' + str(global_sp[line-1][0]) + ' - ' + str(global_sp[line-1][1]) + ')'
            else:
                sp_str = '                     '
            
            if l4 >= line:
                ir_str = ir_str + str(ir_clstrs[line-1][0]) + ' - ' + str(ir_clstrs[line-1][1])
                ir_str = ir_str + ' (' + str(global_ir[line-1][0]) + ' - ' + str(global_ir[line-1][1]) + ')'
            else:
                ir_str = '                    '
            
            if l5 >= line:
                st_str = st_str + str(s_types[line-1][0])
            else:
                st_str = '             '


            line_str = rf_str + pd_str + st_str + sp_str + ir_str
            print(line_str)
        
    
    cur.close()
    db.close()
            
    
    
def firing_order(pri_elems):    #[ [cid cmin cmax], [cid cmin cmax]]
    unique_elems = []
    unique_clstrs =[]
    for row in pri_elems:
        unique_elems = add_unique(row[1:], unique_elems)
    unique_elems.sort()
    #print("unique_elems",unique_elems)
    seq = []
    for row in pri_elems:
        c_num = 1
        for uc in unique_elems:
            if row[1] == uc[0][0]:
                seq.append(c_num)
                break
            c_num +=1
    for row in unique_elems:
        unique_clstrs.append(row[0])
    start_pos = 0
    for elem in seq:
        if elem == 1:
            break
        start_pos += 1
    if start_pos != 0:
        seq = seq[start_pos:] + seq[:start_pos]
    return seq,unique_clstrs

def mt_str(mod_type):
    if mod_type == 'D':
        out_str = 'CONSTANT'
    elif mod_type == 'M':
        out_str = 'STAGGER'
    elif mod_type == 'S':
        out_str = "DWELL/SWITCH"
    elif mod_type == 'J':
        out_str = 'JITTER'
    return out_str
    
    
    
    
    
    
    
    
    
