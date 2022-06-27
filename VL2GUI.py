#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 09:18:05 2022

@author: russelmiller
"""

from tkinter import ttk  #for TreeView
from tkinter import messagebox as msg
import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)


from VL2db_utils import db_con,get_tuning,delete_env
from VL2config import min_baseline_ints,srcdb,fq_refdb
from VL2shared_utils import trunc2,rip
from VL2bb import build_baseline,baseline_parm
from VL2pretty import firing_order
from VL2associator import new_associator
from VL2residue_miner import new_residue_miner
from VL2geo_associator import associate_geo
from VL2geominer import proc_geo
from VL2export import export_pe, make_refdb
from VL2test_utils import int_bldr

import make_VL2pgdb

import configparser
config=configparser.ConfigParser()
config.read('db.config')
refdb_name = config['refdb']['database']

##############################################################################
######                           SHOW/HIDE FCTNS                        ######
##############################################################################

    
def pack_import_wdgts():
    db_import_btn.grid(row=5,column=6)
    db_lagdays_label.grid(row=2,column=4)
    db_lagdays_entry.grid(row=2,column=5)
    db_count_btn.grid(row=4,column=6)

def forget_import_wdgts():
    db_import_btn.grid_forget()
    db_lagdays_label.grid_forget()
    db_lagdays_entry.grid_forget()
    db_count_btn.grid_forget()

def pack_proc_wdgts():
    db_process_btn.grid(row=5,column=6)

def forget_proc_wdgts():
    db_process_btn.grid_forget()
    
def show_db_import():
    hide_db_pnls()
    forget_proc_wdgts()
    pack_import_wdgts()
    dbq_pnl.pack(side=tk.BOTTOM,expand=True)
    
def show_db_proc():
    hide_db_pnls()
    forget_import_wdgts()
    pack_proc_wdgts()
    dbq_pnl.pack(side=tk.TOP,expand=True)
    
def hide_db_pnls():
    dbq_pnl.pack_forget()
    tuning_pnl.pack_forget()
    rand_int_pnl.pack_forget()
    int_hist_pnl.pack_forget()

def show_tuning_pnl():
    # packed during panel definition
    hide_db_pnls()
    tuning_pnl.pack(side=tk.BOTTOM,expand=True)


def show_rand_int_pnl():
    hide_db_pnls()
    rand_int_pnl.pack(side=tk.BOTTOM,expand=True)

def show_int_histo_pnl():
    hide_db_pnls()
    int_hist_pnl.pack(side=tk.BOTTOM)

def show_clstrs():
    clear_histo()   # ensuure obselete histo removed when clstrs are refreshed
    GetGlobalClstrs()    
    clstr_tree.pack()
    clstr_rfrsh_btn.pack()
    clstr_exit_btn.pack()

    
    #clstr_pnl.pack(side=tk.TOP,pady=20)
    clstr_pnl.grid(row=0,pady=15)

def show_new_clstrs():
    # tree already packed with new results (proc_new_dat)
    clear_histo()   # ensuure obselete histo removed when clstrs are refreshed
    clstr_tree.pack()
    #clstr_pnl.pack(side=tk.TOP,pady=20)
    clstr_pnl.grid(row=0,pady=15)

def hide_clstrs():
    #clstr_pnl.pack_forget()
    clstr_pnl.grid_forget()

def show_modes():
    GetModes()
    mode_tree.pack()
    #grps_pnl.pack(side=tk.TOP,pady=20)
    grps_pnl.grid(row=1,pady=15)

def show_new_modes():
    # tree already packed with new results (proc_new_dat)
    mode_tree.pack()
    #grps_pnl.pack(side=tk.TOP,pady=20)
    grps_pnl.grid(row=1,pady=15)
    
def hide_modes():
    hide_mode_parms()
    #grps_pnl.pack_forget()
    grps_pnl.grid_forget()
    
def show_mode_parms():
    #clear_histo()
    #clear_int_histo
    hide_histo()
    parms_exit_btn.pack()
    pri_tree.pack()
    rf_tree.pack()
    pd_tree.pack()
    sp_tree.pack()
    ir_tree.pack()
    parms_pnl.pack()

def hide_mode_parms():
    parms_pnl.pack_forget()

def hide_histo():
    histo_win.pack_forget()
    
def show_geo():
    GetSites()
    site_tree.pack()
    #geo_pnl.pack(pady=20)
    geo_pnl.grid(row=2,pady=15)

def show_new_geo():
    # tree already packed with new results (proc_new_dat)
    site_tree.pack()
    #geo_pnl.pack(pady=20)
    geo_pnl.grid(row=2,pady=15)

def hide_geo():
    #geo_pnl.pack_forget()
    geo_pnl.grid_forget()

def show_histo():
    hide_mode_parms()
    h_subframe.pack(side=tk.TOP)
    histo_win.pack()
    

def clear_histo():
    for widget in h_plot_frame.winfo_children():
        widget.destroy()
    h_subframe.pack_forget()

    
def show_int_histo():
    hide_mode_parms()   
    ih_subframe.pack(side=tk.BOTTOM,pady=20)
    histo_win.pack()
    
    
def clear_int_histo():
    for widget in ih_plot_frame.winfo_children():
        widget.destroy()
    ih_subframe.pack_forget()

    
###############################################################################


###############################################################################
######                        PROCESS FUNCTIONS                          ######
###############################################################################
def clear_tree(tree):
    for i in tree.get_children():
            tree.delete(i)
            
def clear_db_criteria():
    #Clear selections
    db_numdays_entry.delete(0,'end')
    db_lagdays_entry.delete(0,'end')
    db_val1.delete(0,'end')
    db_val2.delete(0,'end')
    db_val3.delete(0,'end')
    db_val4.delete(0,'end')
    db_val5.delete(0,'end')
    field_var1.set(field_options[0])
    field_var2.set(field_options[0])
    field_var3.set(field_options[0])
    field_var4.set(field_options[0])
    field_var5.set(field_options[0])
    criteria_var1.set(field_options[0])
    criteria_var2.set(field_options[0])
    criteria_var3.set(field_options[0])
    criteria_var4.set(field_options[0])
    criteria_var5.set(field_options[0])
    
    min_lat.delete(0,'end')
    max_lat.delete(0,'end')
    min_lon.delete(0,'end')
    max_lon.delete(0,'end')

def GetGlobalClstrs():
    elnot = elnot_var.get()
    mt = mt_var.get()
    parm = parm_var.get()
    constr = ''
    if elnot != 'ALL':
        constr = constr+ " elnot = '" +elnot+ "'"
    if mt !='ALL':
        if constr == '':
            constr = " mod_type ='" +mt+"'"
        else:
            constr = constr + " and mod_type ='" +mt+"'"
    if parm !='ALL':
        if constr == '':
            constr = " parm ='" +parm+"'"
        else:
            constr = constr + " and parm ='" +parm+"'"
    if constr != '':
        constr =' where' + constr
    clear_tree(clstr_tree)
    #histo_subframe.pack()
    
    db = db_con()
    cur = db.cursor()
    cur.execute("SELECT clstr_id,elnot,parm,mod_type,parm_min,parm_max from global_parm_clstrs "+constr+ "order by parm,mod_type,parm_min")
    #if elnot !='ALL':
    #    cur.execute("SELECT clstr_id,elnot,parm,mod_type,parm_min,parm_max from global_parm_clstrs where elnot=%s order by parm,mod_type,parm_min",(elnot,))
    #else:
    #    cur.execute("SELECT clstr_id,elnot,parm,mod_type,parm_min,parm_max from global_parm_clstrs order by elnot,parm,mod_type,parm_min")
    rows = cur.fetchall()  
    print(rows)
    cur.close()
    db.close()
    for row in rows:
        print(row)
        clstr_tree.insert("", tk.END, values=row)

def clstr_selected(event):
    clstr_list = []
    elnot=None
    parm=None
    mod_type=None
    min_val = 1e6
    max_val = -1
    first_pass = True
    
    for selected_item in clstr_tree.selection():
        item = clstr_tree.item(selected_item)
        record = item['values']
        if first_pass:
            elnot = record[1]
            parm = record[2]
            mod_type = record[3]
            first_pass = False
        if record[1] != elnot or record[2] != parm or record[3] != mod_type:
            msg.showinfo("abort","ELNOT, PARM and MOD_TYPE must be the same for multiple selections!")
            return
        cmin = float(record[4])
        cmax = float(record[5])
        clstr_list.append([cmin,cmax])
        if cmin < min_val:
            min_val = cmin
        if cmax > max_val:
            max_val = cmax

    if elnot != None and parm != None:
        info = [elnot,parm,mod_type,min_val,max_val]
        plot_clstr_histo(info,clstr_list)

def plot_clstr_histo(info,clstr_list):
    clear_histo()
    show_histo()
    
    elnot=info[0]
    parm=info[1]
    mod_type = info[2]
    pmin=float(info[3])
    pmax=float(info[4])
    #print('elnot',elnot)
    #print('mod_type',mod_type)
    #print('parm',parm)
    #print('seperator:','/')
    p_title=elnot+'/'+mod_type+'/'+parm
    
    #### Buffer by +/- 2 horizons
    incr,horizon,p = get_tuning(elnot,parm)
    parm_min = max(0,pmin-2*horizon)
    parm_max = pmax+2*horizon

    #### Retrieve intercept data
    db = db_con()
    cur = db.cursor()
    
    if parm.upper() =='PRI':
        q_str = """select pri_value from intercepts a,intercept_pris b where a.intercept_id = b.intercept_id and elnot=%s 
                    and adj_mt(a.mod_type)=%s and pri_value between %s and %s"""
    else:
        q_str = "select " +parm+ " from intercepts where elnot = %s and adj_mt(mod_type)=%s and " +parm+ " between %s and %s"
    cur.execute(q_str,(elnot,mod_type,parm_min,parm_max)) 
    dat = cur.fetchall()
    
    if dat == []:
        print("No "+parm+" intercept data exists in the range "+str(parm_min)+" - "+str(parm_max))
        return
    else:
        print(len(dat),"values retrieved")
    
    dat = np.array(dat).astype(float)
    
    n_bins = min(round((parm_max-parm_min)/incr),1000)
    print("nbins",n_bins)
    
    # the figure that will contain the plot
    fig = Figure(figsize = (5, 5),dpi = 100)
  
     
    # adding the plots to the figure
    plot1 = fig.add_subplot(111,title=p_title)
    hist = plot1.hist(dat,bins = n_bins)  # plot the histogram
    y=max(hist[0])
    for row in clstr_list:
        pmin = row[0]
        pmax = row[1]
        plot1.plot([pmin,pmin,pmax,pmax],[0,y,y,0],color='red')  # overlay the clustr boundaries
  

    # creating the Tkinter canvas
    
    canvas = FigureCanvasTkAgg(fig,master = h_plot_frame)
    #canvas = FigureCanvasTkAgg(fig,master = histo_win)
    canvas.draw()
  
    # placing the canvas on the Tkinter window
    canvas.get_tk_widget().pack(side=tk.LEFT)

def pop_histo():
    x_win=tk.Tk()
    x_frame = tk.Frame(x_win)
    x_frame.pack()
    # Grab selected data & redraw in external window - mostly copy/paste from clstr_selected and plot_clstr_histo
    clstr_list = []
    min_val = 1e6
    max_val = -1
    first_pass = True
    for selected_item in clstr_tree.selection():
        item = clstr_tree.item(selected_item)
        record = item['values']
        if first_pass:
            elnot = record[1]
            parm = record[2]
            mod_type = record[3]
            first_pass = False
        if record[1] != elnot or record[2] != parm or record[3] != mod_type:
            msg.showinfo("abort","ELNOT, PARM and MOD_TYPE must be the same for multiple selections!")
            return
        cmin = float(record[4])
        cmax = float(record[5])
        clstr_list.append([cmin,cmax])
        if cmin < min_val:
            min_val = cmin
        if cmax > max_val:
            max_val = cmax
    info = [elnot,parm,mod_type,min_val,max_val]
    elnot=info[0]
    parm=info[1]
    mod_type = info[2]
    pmin=float(info[3])
    pmax=float(info[4])
    p_title=elnot+'/'+mod_type+'/'+parm
    
    #### Buffer by +/- 2 horizons
    incr,horizon,p = get_tuning(elnot,parm)
    parm_min = max(0,pmin-2*horizon)
    parm_max = pmax+2*horizon

    #### Retrieve intercept data
    db = db_con()
    cur = db.cursor()
    
    if parm.upper() =='PRI':
        q_str = """select pri_value from intercepts a,intercept_pris b where a.intercept_id = b.intercept_id and elnot=%s 
                    and adj_mt(a.mod_type)=%s and pri_value between %s and %s"""
    else:
        q_str = "select " +parm+ " from intercepts where elnot = %s and adj_mt(mod_type)=%s and " +parm+ " between %s and %s"
    cur.execute(q_str,(elnot,mod_type,parm_min,parm_max)) 
    dat = cur.fetchall()
    
    if dat == []:
        print("No "+parm+" intercept data exists in the range "+str(parm_min)+" - "+str(parm_max))
        return
    else:
        print(len(dat),"values retrieved")
    
    dat = np.array(dat).astype(float)
    
    n_bins = min(round((parm_max-parm_min)/incr),1000)
    print("nbins",n_bins)
    
    # the figure that will contain the plot
    fig = Figure(figsize = (5, 5),dpi = 100)
  
     
    # adding the plots to the figure
    plot1 = fig.add_subplot(111,title=p_title)
    hist = plot1.hist(dat,bins = n_bins)  # plot the histogram
    y=max(hist[0])
    for row in clstr_list:
        pmin = row[0]
        pmax = row[1]
        plot1.plot([pmin,pmin,pmax,pmax],[0,y,y,0],color='red')  # overlay the clustr boundaries
  

    # creating the Tkinter canvas
    
    canvas = FigureCanvasTkAgg(fig,master = x_frame)
    #canvas = FigureCanvasTkAgg(fig,master = histo_win)
    canvas.draw()
  
    # placing the canvas on the Tkinter window
    canvas.get_tk_widget().pack()
    x_win.mainloop()

def plot_int_histo():
    clear_int_histo()
    show_int_histo()
    elnot = elnot_var.get()
    parm = parm_var.get()
    mt = mt_var.get()
    if elnot =='ALL' or parm =='ALL':
        msg.showinfo("warn","Please select desired parameters from the pull-downs")
        return
    p_min = ih_min_entry.get()
    p_max = ih_max_entry.get()
    try:
        p_min = float(p_min)
    except:
        msg.showinfo("warn", "Min val must be a real number")
        return
    try:
        p_max = float(p_max)
    except:
        msg.showinfo("warn", "Max val must be a real number")
        return
    int_histo(elnot,parm,mt,p_min,p_max)
    #int_histo_subframe.pack()
    
    
def int_histo(elnot,parm,mod_type,parm_min,parm_max):
    incr,horizon,p = get_tuning(elnot,parm)  # incr for number of bins calc
    #pmin=parm_min #end points for red overlay Note: since no longer adding 2H buffer, overlay not really needed...
    #pmax=parm_max
    
    p_title=elnot+'/'+mod_type+'/'+parm
    
    if parm.upper() == 'PRI':
        dbq = "select pri_value from intercepts a,intercept_pris b where a.intercept_id=b.intercept_id and elnot=%s "
    else:
        dbq = "select " +parm.lower()+ " from intercepts where elnot=%s "
        
    if mod_type !='ALL':
        dbq = dbq + "and adj_mt(mod_type)=%s "
    
    if parm.upper() == 'PRI':
        dbq = dbq + "and pri_value between %s and %s"
    else:
        dbq = dbq + "and " +parm.lower()+ " between %s and %s"
    
    db = db_con()
    cur = db.cursor()
    
    if mod_type == 'ALL':
        cur.execute(dbq,(elnot,parm_min,parm_max))
    else:
        cur.execute(dbq,(elnot,mod_type,parm_min,parm_max))   
    
    dat = cur.fetchall()
    if dat == []:
        err_msg = "No "+parm+" intercept data exists in the range "+str(parm_min)+" - "+str(parm_max)
        msg.showinfo("warn",err_msg)
        return
    else:
        print(len(dat),"values retrieved")
    
    dat = np.array(dat).astype(float)
    
    n_bins = min(round((parm_max-parm_min)/incr),1000)
    print("nbins",n_bins)
    
    # the figure that will contain the plot
    fig = Figure(figsize = (5, 5),dpi = 100)
  
     
    # adding the plots to the figure
    plot1 = fig.add_subplot(111,title=p_title)
    hist = plot1.hist(dat,bins = n_bins)  # plot the histogram
    
    # since histo limits not buffered by 2 x horizon, the limits overlay is counterproductive
    #y=max(hist[0])
    #plot1.plot([pmin,pmin,pmax,pmax],[0,y,y,0],color='red')  # overlay the clustr boundaries
  
    # creating the Tkinter canvas
    
    canvas = FigureCanvasTkAgg(fig,master = ih_plot_frame)  
    canvas.draw()
  
    # placing the canvas on the Tkinter window
    canvas.get_tk_widget().pack()

def GetModes():
    elnot = elnot_var.get()
    mt = mt_var.get()
    constr = ''
    if elnot != 'ALL':
        constr = constr+ " elnot = '" +elnot+ "'"
    if mt !='ALL':
        if constr == '':
            constr = " mod_type ='" +mt+"'"
        else:
            constr = constr + " and mod_type ='" +mt+"'"
    
    if constr != '':
        constr =' where' + constr
    clear_tree(mode_tree)
    clear_tree(pri_tree)
    clear_tree(rf_tree)
    clear_tree(pd_tree)
    clear_tree(sp_tree)
    clear_tree(ir_tree)

    db = db_con()
    cur = db.cursor()
    cur.execute("SELECT group_id,elnot,mod_type,group_descr from parameter_groups "+constr+" order by elnot,mod_type")
    
    #if elnot !='ALL':
    #    cur.execute("SELECT group_id,elnot,mod_type,group_descr from parameter_groups where elnot=%s order by mod_type",(elnot,))
    #else:
    #    cur.execute("SELECT group_id,elnot,mod_type,group_descr from parameter_groups order by elnot,mod_type")
    dat = cur.fetchall()
    rows = []
    for row in dat:
        group_id = row[0]
        elnot = row[1]
        mod_type = row[2]
        group_descr = row[3]
        
        cur.execute(""" select a.clstr_id, b.parm_min,b.parm_max from sequence_elements a, global_parm_clstrs b, sequences c 
                    where a.clstr_id = b.clstr_id  and a.seq_id =c.seq_id and c.group_id = %s and c.is_parent = 'T' order by a.position""",(group_id,))
        pri_elems = cur.fetchall()
        seq,seq_clstrs = firing_order(pri_elems)
        seq_str = ''
        for elem in seq:
            if seq_str != '':
                seq_str = seq_str + '-'
            seq_str = seq_str + str(elem)
        rows.append([group_id,elnot,mod_type,seq_str,group_descr])
    cur.close()
    db.close()
    for row in rows:
        mode_tree.insert("", tk.END, values=row)

def mode_selected(event):
    group_id = 0
    for selected_item in mode_tree.selection():
        item = mode_tree.item(selected_item)
        record = item['values']
        group_id = int(record[0])
        GetModeParms(group_id)

def GetModeParms(group_id):
    clear_tree(pri_tree)
    clear_tree(rf_tree)
    clear_tree(pd_tree)
    clear_tree(sp_tree)
    clear_tree(ir_tree)
    
    #l = len(pri_tree.get_children())
    #if l>0:
    #    for i in pri_tree.get_children():
    #        pri_tree.delete(i)
    
    db = db_con()
    cur = db.cursor()
    print("group_id =",group_id)
    cur.execute("""select distinct a.parm_min,a.parm_max from global_parm_clstrs a, sequences b, sequence_elements c
                         where a.parm='PRI' and a.clstr_id=c.clstr_id and b.seq_id=c.seq_id and b.group_id = %s order by a.parm_min""",
                         (group_id,))
    rows = cur.fetchall()
    print(rows)
    i = 1
    for row in rows:
        row2 = (i,) + row
        i+=1
        print(row2)
        pri_tree.insert("", tk.END, values=row2)
    
    #RF
    l = len(rf_tree.get_children())
    if l>0:
        for i in rf_tree.get_children():
            rf_tree.delete(i)
    cur .execute("""select a.clstr_id,a.parm_min,a.parm_max,b.parm_min,b.parm_max from global_parm_clstrs a, local_parm_clstrs b where
                 a.parm = 'RF' and a.parm=b.parm and a.clstr_id=b.clstr_id and b.group_id = %s order by a.parm_min""",(group_id,))
    rows = cur.fetchall()
    for row in rows:
        rf_tree.insert("", tk.END, values=row)
    
    #PD
    l = len(pd_tree.get_children())
    if l>0:
        for i in pd_tree.get_children():
            pd_tree.delete(i)
    cur .execute("""select a.clstr_id,a.parm_min,a.parm_max,b.parm_min,b.parm_max from global_parm_clstrs a, local_parm_clstrs b where
                 a.parm = 'PD' and a.parm=b.parm and a.clstr_id=b.clstr_id and b.group_id = %s order by a.parm_min""",(group_id,))
    rows = cur.fetchall()
    for row in rows:
        pd_tree.insert("", tk.END, values=row)
    
    #SP
    l = len(sp_tree.get_children())
    if l>0:
        for i in sp_tree.get_children():
            sp_tree.delete(i)
    cur .execute("""select a.clstr_id,a.parm_min,a.parm_max,b.parm_min,b.parm_max from global_parm_clstrs a, local_parm_clstrs b where
                 a.parm = 'SP' and a.parm=b.parm and a.clstr_id=b.clstr_id and b.group_id = %s order by a.parm_min""",(group_id,))
    rows = cur.fetchall()
    for row in rows:
        sp_tree.insert("", tk.END, values=row)
    
    #IR
    l = len(ir_tree.get_children())
    if l>0:
        for i in ir_tree.get_children():
            ir_tree.delete(i)
    cur .execute("""select a.clstr_id,a.parm_min,a.parm_max,b.parm_min,b.parm_max from global_parm_clstrs a, local_parm_clstrs b where
                 a.parm = 'IR' and a.parm=b.parm and a.clstr_id=b.clstr_id and b.group_id = %s order by a.parm_min""",(group_id,))
    rows = cur.fetchall()
    for row in rows:
        ir_tree.insert("", tk.END, values=row)
    
    cur.close()
    db.close()
    show_mode_parms()

def GetSites():
   
    clear_tree(site_tree)
    elnot = elnot_var.get()
    db = db_con()
    cur = db.cursor()
    if elnot != 'ALL':
        cur.execute("SELECT site_id,elnot,sma,smi,orient,lat,lon,heard_count,roa from site_geos where elnot=%s",(elnot,))
    else:
        cur.execute("SELECT site_id,elnot,sma,smi,orient,lat,lon,heard_count,roa from site_geos order by elnot")
    rows = cur.fetchall()
    cur.close()
    db.close()

    for row in rows:
        site_tree.insert("", tk.END, values=row)
    
def count_vipr_ints():
    db = db_con()
    cur = db.cursor()
    clause = build_db_clause(True)
    
    q_str = "select count(*) from " + srcdb + ".intercepts"
    if clause != '':
        q_str = q_str + " where " + clause
    print(q_str)
    cur.execute(q_str)
    cnt = cur.fetchone()[0]
    msg.showinfo("info","Intercept query returns " + str(cnt) + " rows")
    cur.close()
    db.close()

def build_db_clause(t0_is_now):
    start_delta = None
    stop_delta = None
    tp_clause = ''
    p_elnot = elnot_var.get()
    num_days = db_numdays_entry.get()
    valid,vmsg = validate(num_days,'I',3)
    if not valid:
        num_days = ''
        msg.showwarning(title=None, message='num_days input ignored: '+vmsg)
    lag_days = db_lagdays_entry.get()
    valid,vmsg = validate(lag_days,'I',3)
    if not valid:
        lag_days = ''
        msg.showwarning(title=None, message='lag_days input ignored: '+vmsg)
    
    #### Pull-Down Constraints
    if p_elnot not in ('','ALL'):
        clause = "elnot='" +p_elnot + "'"
    else:
        clause = ''            
    for i in range(1,6):        
        f = eval("field_var" + str(i) +".get()")
        c = eval("criteria_var" + str(i) +".get()")
        v = eval("db_val" + str(i) +".get()")
        
        # Validate string entries retrieved into v
        if v != '' and f != '':
            if f.upper() == 'ELNOT':
                valid,vmsg = validate(v,'AN',5)
            elif f.upper() in ('MOD_TYPE','COUNTRY_CODE','COLLECTOR','RD_OUT_STAT'):
                valid,vmsg = validate(v,'A',2)
            elif f.upper() == 'NUM_BURSTS':
                valid,vmsg = validate(v,'I',10)
            else:
                msg.showwarning(title=None, message=f+' has not been validated...contact developer')
            if not valid:
                v = ''
                msg.showwarning(title=None, message=f+' input ignored: '+vmsg)
            
        
        # end validation
        if f != '' and c != '' and v != '':
            if clause != '':
                clause = clause + " and "
            if c in ('IN','NOT IN'):
                clause = clause + f + " " + c + list_clause(v.upper())
            else:
                clause = clause + f + c + "'" + v.upper() + "'"
            #print('intermediate clause',clause)     
    #### Find geo_constraints
    lat_min = min_lat.get()
    valid,vmsg = validate(lat_min,'D',10)
    if not valid:
        lat_min = ''
        msg.showwarning(title=None, message='lat_min input ignored: '+vmsg)
    lat_max = max_lat.get()
    valid,vmsg = validate(lat_max,'D',10)
    if not valid:
        lat_max = ''
        msg.showwarning(title=None, message='lat_max input ignored: '+vmsg)
    lon_min = min_lon.get()
    valid,vmsg = validate(lon_min,'D',10)
    if not valid:
        lon_min = ''
        msg.showwarning(title=None, message='lon_min input ignored: '+vmsg)
    lon_max = max_lon.get()
    valid,vmsg = validate(lon_max,'D',10)
    if not valid:
        lon_max = ''
        msg.showwarning(title=None, message='lon_max input ignored: '+vmsg)
    
    geo_list = []
    if lat_min != '':
        geo_list.append(" latitude > " + lat_min)
    if lat_max != '':
        geo_list.append(" latitude < " + lat_max)
    
    if lon_min != '':
        geo_list.append(" longitude > " + lon_min)
    if lon_max != '':
        geo_list.append(" longitude < " + lon_max)
        
    geo_clause = ''
    for sub_clause in geo_list:
        if geo_clause == '':
            geo_clause = sub_clause
        else:
            geo_clause = geo_clause + " and " + sub_clause
    
    if geo_clause != '':
        if clause == '':
            clause = geo_clause
        else:
            clause = clause + " and " + geo_clause
        
    

                    
    #### find time_process constraints
    if num_days != '':
        if lag_days != '':
            start_delta = num_days
            stop_delta = lag_days        
        else:
            start_delta = num_days            
    if t0_is_now:
        if start_delta != None:
            tp_clause = "time_process > current_timestamp - interval '" +start_delta+ "' day"
        if stop_delta != None:
            tp_clause = tp_clause + " and time_process < current_timestamp -interval '" +stop_delta+ "' day"
        if clause != '' and tp_clause != '':
            clause = clause + " and "
        clause = clause + tp_clause
        #print('tp_clause',tp_clause)
        return clause    
    else:  # starting from a previous time
        db = db_con()
        cur = db.cursor()
        max_id_str = "select max(intercept_id) from intercepts"
        if clause != '':
            max_id_str = max_id_str + " where " + clause 
        cur.execute(max_id_str)
        max_id = cur.fetchone()[0]
        if max_id != None:
            cur.execute("select time_process from intercepts where intercept_id=" + str(max_id))
            start_proc = cur.fetchone()[0]
        
            tp_clause = " intercept_id > " +str(max_id) + " and time_process < timestamp '" +str(start_proc) + "' + interval '" + num_days + "' day"
        #print('pre-clause',clause)
        #print("tp_clause",tp_clause)
        if clause != '' and tp_clause != '':
            clause = clause + " and " + tp_clause
        elif tp_clause != '':
            clause = tp_clause
        cur.close()
        db.close()
        return clause
    return ''

def refresh_elnot_list():
    db = db_con()
    cur= db.cursor()
    cur. execute("select distinct elnot from intercepts order by elnot")
    dat = cur.fetchall()
    elnot_options=['ALL']
    for notation in dat:
        elnot_options.append(notation[0])
    cur.close()
    db.close()
    menu = elnot_pull_down["menu"]
    menu.delete(0, "end")
    for string in elnot_options:
        menu.add_command(label=string, 
                             command=lambda value=string: elnot_var.set(value))
    if len(elnot_options)>1:
        elnot_var.set(elnot_options[1])
    else:
        elnot_var.set(elnot_options[0])

def list_clause(v):
    out_str = " ('"
    for char in v:
        if char == ',':
            out_str = out_str + "','"
        else:
            if char != ' ':   #ignore blank spaces
                out_str = out_str + char
    out_str = out_str + "')"
    return out_str


def validate(list_value,item_type,max_item_len):
    # Validate a user string entry of an item or csv list
    # item_type: 'A' (alpha), 'AN' (alphanumeric), 'D' (decimal value), or 'I' (integer value)
    # max_item_len: maximum length of each individual item
    
    list_item = ''
    for char in list_value:
        if char == ',':
            valid,reason = validate_item(list_item,item_type,max_item_len) #validate the list item preceeding the comma
            if valid:
                list_item = ''
            else:
                return valid,reason
        else:
            if char != ' ':     #ignore blank spaces
                list_item = list_item + char
    
    valid,reason = validate_item(list_item,item_type,max_item_len) #validate the last (or only) item in the list
    return valid,reason

def validate_item(value,value_type,max_len):
    # value: a string (entered by a user) to be validated
    # value_type: A (Alpha only), AN = (AlphaNumeric), D (decimal) or I (integer)
    # max_len: the maximum length of the input string
    #
    # Returns True and a 'success' string if value is validated
    # Returns False and an 'error' string if value fails validation
    if value == '':
        return True, "empty string"
    if len(value) > max_len:
        return False, 'Value provided contains more characters than alowed'
    if value_type.upper() == 'A':
        for i in value:
            if i.isalpha():
                continue
            return False,"Input contains non alphabetic characters"
    elif value_type.upper() == 'AN':
        for i in value:
            if i.isalnum():
                continue
            return False,"Input contains non-alphanumeric characters"
    elif value_type.upper() == 'D':
        try:
            test_val = float(value)
        except:
            return False, "Only decimal values allowed"
    elif value_type == 'I':
        try:
            test_val = int(value)
        except:
            return False, "Only integer values allowed"
    return True, "Input validated"

def import_vipr_ints():
    db = db_con()
    cur = db.cursor()
    cur.execute("select count(*) from intercepts")
    int_cnt = cur.fetchone()[0]
    cur.execute("select count(*) from intercept_pris")
    pri_cnt = cur.fetchone()[0]
    
    clause = build_db_clause(True)
    q_str = "insert into intercepts select * from " +srcdb+ ".intercepts "
    if clause != '':
        q_str = q_str + " where " + clause + " and intercept_id not in (select intercept_id from intercepts)"
    else:
        q_str = q_str + "where intercept_id not in (select intercept_id from intercepts)"
    cur.execute(q_str)
    print(q_str)
    
    ####
    q_str = """insert into intercept_pris select * from """ +srcdb+ """.intercept_pris where intercept_id in (select intercept_id from intercepts)
               and intercept_id not in (select intercept_id from intercept_pris)"""
    cur.execute(q_str)
    
    print(q_str)
    
    cur.execute("commit")
    cur.execute("select count(*) from intercepts")
    int_cnt2 = cur.fetchone()[0]
    cur.execute("select count(*) from intercept_pris")
    pri_cnt2 = cur.fetchone()[0]
    
    msg_str = str(int_cnt2-int_cnt) + " intercepts loaded with " + str(pri_cnt2-pri_cnt) + " PRIs"
    print(msg_str)
    
    # Refresh elnot_options:
    refresh_elnot_list()
    msg.showinfo("result",msg_str)
    
    cur.close()
    db.close()
    
def delete_VL2_ints():
    x = msg.askokcancel("warn","Delete all VL2 intercepts and the entire preserved environment?")
    if not x:
        return
    db = db_con()
    cur = db.cursor()
    cur.execute("delete from intercept_pris")
    cur.execute("delete from intercepts")    
    cur.execute("commit")
    
    # delete the PE
    delete_env(0)
    # Refresh elnot_options:
    refresh_elnot_list()
        
    msg.showinfo("done","All VL2 intercepts and the PE deleted!")
    cur.close()
    db.close()
    
def delete_VL2_elnot_ints():
    d_elnot = elnot_var.get()
    if d_elnot =='ALL':
        msg.showinfo("error","Error - No ELNOT selected")
        return
    x = msg.askokcancel("warn","Delete all "+d_elnot+" intercepts and corresponding preserved environment?")
    if not x:
        return
    db = db_con()
    cur = db.cursor()
    cur.execute("delete from intercept_pris where intercept_id in (select intercept_id from intercepts where elnot=%s)",(d_elnot,))
    cur.execute("delete from intercepts where elnot=%s",(d_elnot,))    
    cur.execute("commit")
    cur.close()
    db.close()
    refresh_elnot_list()
    delete_env(d_elnot)
    msg_str = d_elnot + "intercepts and associated PE deleted!"
    msg.showinfo("done",msg_str)

def update_VL2_schema():
    make_VL2pgdb
    msg.showinfo("done","VL2 DB Schema Updated!")

def QuitApp():
    root.destroy()

def rebaseline():
    # Rebuild parametric baseline after cluster tuning
    elnot = elnot_var.get()
    if elnot =='ALL':
        msg.showinfo("error","Error - No ELNOT selected")
        return
    msg_str="Re-baseline " +elnot+ "?"
    x = msg.askokcancel("Are you sure?",msg_str)
    # alert :Reprocess BBBBB?
    if x:
        build_baseline(elnot,False)
    msg_str = "Parametric Baseline for " +elnot+ " is complete"
    msg.showinfo("done",msg_str)
    show_clstrs()
    show_modes()
    show_geo()

def re_clstr():
    print(1)
    
    elnot = elnot_var.get()
    if elnot =='ALL':
        msg.showinfo("error","Error - No ELNOT selected")
        return
    parm = parm_var.get()
    ok=msg.askokcancel("sure","Reprocess "+elnot+ " " +parm+"?")
    if not ok:
        return
    delete_env(elnot,'P',False)  # only parametrics, don't delete global clstrs or associated residue
    print(3)
    cases = delete_global_parm(elnot,parm)
    for mod_type in cases:
        print(mod_type)
        baseline_parm(elnot,mod_type,parm)
    msg.showinfo("done","re-clustering complete")
    show_clstrs()
    show_modes()

def delete_global_parm(elnot,parm):
    db = db_con()
    cur = db.cursor()
    cur.execute("select distinct mod_type from global_parm_clstrs where elnot=%s and parm=%s",(elnot,parm))
    cases = rip(cur.fetchall(),0)
    cur.execute("delete from global_parm_clstrs where elnot=%s and parm=%s",(elnot,parm))
    cur.execute("delete from parm_residue where elnot=%s and parm=%s",(elnot,parm))
    cur.execute("commit")
    
    cur.close()
    db.close()
    return cases

def gen_elnot_baseline():
    b_elnot=elnot_var.get()
    if b_elnot =='ALL':
        msg.showinfo("error","Error - No ELNOT selected")
        return
    build_baseline(b_elnot)
    msg_str = "Baseline for " +b_elnot+ " Complete!"
    msg.showinfo("done",msg_str)
    show_clstrs()
    show_modes()
    show_geo()

def gen_full_baseline():
    elnots = get_int_elnots()
    if elnots == []:
        msg.showinfo("warn","No data to process...load some intercepts and try again")
        return
    for b_elnot in elnots:
        build_baseline(b_elnot)
    msg.showinfo("done","New PE Baseline Complete!")
    show_clstrs()
    show_modes()
    show_geo()

def get_int_elnots():
    db = db_con()
    cur=db.cursor()
    cur.execute("select elnot,count(elnot) from intercepts group by elnot having count(elnot)>%s order by elnot",(min_baseline_ints,))
    elnots = rip(cur.fetchall(),0)
    cur.close()
    db.close()
    return elnots

def delete_local_env():
    x = msg.askokcancel("Confirm", "Are you sure you want to delete the VL2 Preserved Environment?")
    if x:
        delete_env(0)
        msg.showinfo("done",'Local VL2 PE Deleted!')
        show_clstrs()
        show_modes()
        show_geo()
def proc_new_dat():
    #db_lagdays_entry.delete(0,'end')
    num_days = db_numdays_entry.get()
    valid,vmsg = validate(num_days,'I',3)
    if not valid:
        num_days = ''
        msg.showwarning(title=None, message='num_days input ignored: '+vmsg)
        
    p_elnot=elnot_var.get()
    if num_days == '':
        msg.showinfo("error","Number of days must be specified...Update aborted")
        return
    clause = build_db_clause(False)
    print('final clause',clause)
    
    db=db_con()
    cur=db.cursor()
    
    # Ensure data is available
    db_count_query = "select count(*) from " +srcdb+ ".intercepts where " + clause
    cur.execute(db_count_query)
    int_count = cur.fetchone()[0]
    print('int_count',int_count)
    if int_count == 0:
        msg.showinfo("abort","No data available in the specified range: "+clause)
        return
    else:
        ok=msg.askokcancel("verify","Query for: "+clause+" returns "+str(int_count)+" rows...Continue?")
        if not ok:
            return

    # Capture pre-processing stats
    cur.execute("select max(clstr_id) from global_parm_clstrs")
    pre_clstr_max = cur.fetchone()[0]
    
    cur.execute("select max(group_id) from parameter_groups")
    pre_grp_max = cur.fetchone()[0]
    
    cur.execute("select max(site_id) from site_geos")
    pre_site_max = cur.fetchone()[0]
    
    
    # Load & stage Data
    cur.execute("delete from comparison_queue")
    db_stage_query = "insert into comparison_queue(intercept_id) select intercept_id from " +srcdb+ ".intercepts where " + clause
    db_int_query = "insert into intercepts select * from " +srcdb+ ".intercepts where intercept_id in (select intercept_id from comparison_queue)"
    db_int_pri_query = "insert into intercept_pris select * from " +srcdb+ ".intercept_pris where intercept_id in (select intercept_id from comparison_queue)"
    print(db_stage_query)
    print(db_int_query)
    print(db_int_pri_query)
    cur.execute(db_stage_query)
    cur.execute(db_int_query)
    cur.execute(db_int_pri_query)
    cur.execute("commit")
    
    new_associator()
    associate_geo()
    if p_elnot == 'ALL':
        new_residue_miner()
    else:
        new_residue_miner(p_elnot)
    cur.execute("select distinct a.elnot from intercepts a, residue_ints b where a.intercept_id=b.intercept_id and b.r_type= 'G'")
    cases = cur.fetchall()
    for case in cases:
        proc_geo(case)
    
    # Capture new results
    #   - Basic for now, but eventually display new results in tree structures
    cur.execute("select clstr_id,elnot,parm,mod_type,parm_min,parm_max from global_parm_clstrs where clstr_id > %s",(pre_clstr_max,))
    new_clstrs = cur.fetchall()
    new_clstr_cnt = len(new_clstrs)
    print("new_clstrs:",new_clstrs)
    
    cur.execute("select group_id,elnot,mod_type,group_descr from parameter_groups where group_id > %s",(pre_grp_max,))
    dat = cur.fetchall()
    new_grp_cnt = len(dat)
    #print("new_groups:",new_groups)
    
    ##
    new_groups = []
    for row in dat:
        group_id = row[0]
        elnot = row[1]
        mod_type = row[2]
        group_descr = row[3]
        
        cur.execute(""" select a.clstr_id, b.parm_min,b.parm_max from sequence_elements a, global_parm_clstrs b, sequences c 
                    where a.clstr_id = b.clstr_id  and a.seq_id =c.seq_id and c.group_id = %s and c.is_parent = 'T' order by a.position""",(group_id,))
        pri_elems = cur.fetchall()
        seq,seq_clstrs = firing_order(pri_elems)
        seq_str = ''
        for elem in seq:
            if seq_str != '':
                seq_str = seq_str + '-'
            seq_str = seq_str + str(elem)
        new_groups.append([group_id,elnot,mod_type,seq_str,group_descr])
        
    ##
    
    
    cur.execute("select site_id,elnot,sma,smi,orient,lat,lon,heard_count,roa from site_geos where site_id > %s",(pre_site_max,))
    new_sites = cur.fetchall()
    new_site_cnt = len(new_sites)
    print("new_sites:",new_sites)
    
    cur.close()
    db.close()
    result_str = str(new_clstr_cnt) + " clusters; " + str(new_grp_cnt) + "modes; and " + str(new_site_cnt) + " sites found. Processing complate!"
    
    # Populate trees with new data
    clear_tree(clstr_tree)
    clear_tree(mode_tree)
    clear_tree(site_tree)
    clear_tree(pri_tree)
    clear_tree(rf_tree)
    clear_tree(pd_tree)
    clear_tree(sp_tree)
    clear_tree(ir_tree)
    
    
    for row in new_clstrs:
        #print(row)
        clstr_tree.insert("", tk.END, values=row)
    
    for row in new_groups:
        #print(row)
        mode_tree.insert("", tk.END, values=row)
    
    for row in new_sites:
        #print(row)
        site_tree.insert("", tk.END, values=row)
    
    
    msg.showinfo("done",result_str)
    show_new_clstrs()
    show_new_modes()
    show_new_geo()

def get_current_parms():
    elnot = elnot_var.get()
    if elnot =='ALL':
        msg.showinfo("error","Error - No ELNOT selected")
        return
    parm = parm_var.get()  #[2:5]
    if elnot == 'ALL' or parm == 'ALL':
        msg.showinfo("warn","You must select ELNOT & parameter in the pull-down menus")
        return
    incr,horizon,thresh,p = get_tuning(elnot,parm,'BOTH')
    out_str = str(elnot) + " - " +str(parm)+ ":    Incr = " +str(incr)+ ",  Horizon = " +str(horizon)+ ",  Threshold = " +str(thresh)+ ",  Percentile = " +str(p)
    msg.showinfo("Current Tuning Parms",out_str)
    print(out_str)
    
def update_parms():
    x = False
    parm_list = []
    new_h = h_entry.get()
    valid,vmsg =validate(new_h,'D',5)
    if not valid:
        new_h = ''
        msg.showwarning(title='Input Error',message = 'horizon input ignored: '+vmsg)
    new_t = t_entry.get()
    valid,vmsg =validate(new_t,'D',5)
    if not valid:
        new_t = ''
        msg.showwarning(title='Input Error',message = 'threshold input ignored: '+vmsg)
    
    new_p = p_entry.get()
    valid,vmsg =validate(new_p,'D',5)
    if not valid:
        new_p = ''
        msg.showwarning(title='Input Error',message = 'percentile input ignored: '+vmsg)
    
    if new_h == '':
        h_str = ''
    else:
        h_str = 'Horizon = ' +new_h
        parm_list.append(['horizon',new_h])
    if new_t == '':
        t_str = ''
    else:
        t_str = '  Threshold = ' +new_t
        parm_list.append(['thresh',new_t])
    if new_p == '':
        p_str = ''
    else:
        p_str = '  p = ' +new_p
        parm_list.append(['p',new_p])
    out_str = h_str + t_str + p_str
    if out_str == '':
        msg.showinfo("No change","No changes made")
    else:
        elnot = elnot_var.get()
        if elnot =='ALL':
            msg.showinfo("error","Error - No ELNOT selected")
            return
        parm = parm_var.get()
        verify_str = "Commit update("+elnot+", " +parm+ "):" +out_str+"?" 
        x = msg.askokcancel("Change",verify_str)
    if x:
        db = db_con()
        cur = db.cursor()
        cur.execute("select count(1) from direct_tune_parms where elnot=%s and parm=%s",(elnot,parm))
        cnt = cur.fetchone()[0]
        
        if cnt == 0:   #insert row
            cur.execute("insert into direct_tune_parms(elnot,parm) values(%s,%s)",(elnot,parm))
        for k in parm_list:
            update_str = "update direct_tune_parms set " + k[0] + " = " + k[1] + "where elnot = %s and parm = %s"
            print("update_str = ",update_str)
            cur.execute(update_str,(elnot,parm))
        cur.execute("commit")
        cur.close()
        db.close()
        print("proceed with update...")

def build_list(v,v_type):
    # v_type = S(string), N(numeric) or I(integer)
    out_lst = []
    item = ''
    for char in v:
        if char == ',':
            if v_type == 'S':
                out_lst.append(item.upper())
            elif v_type == 'N':
                out_lst.append(float(item))
            elif v_type == 'I':
                out_lst.append(int(item))
            item = ''
        else:
            item = item + char
    if v_type == 'S':
        out_lst.append(item.upper())   # capture last item in the list
    if v_type == 'N':
        out_lst.append(float(item))   # capture last item in the list
    if v_type == 'I':
        out_lst.append(int(item))   # capture last item in the list
    return out_lst


def gen_rand_ints():
    delete_var = not append_var.get()
    seeds = []
    try:
        num_ints = int(num_rand_ints.get())
    except:
        msg.showinfo("warn","A valid (integer) number of intercepts must be specified")
        return
    
    #### Retrieve/validate seeds
    
    elnots = elnot_seeds.get()  
    if len(elnots) == 0:
        msg.showinfo("warn","At least one ELNOT seed must be provided")
        return
    valid,vmsg = validate(elnots,'AN',5)
    if not valid:
        msg.showinfo("warn",'ELNOT: '+vmsg)
        return
    
    mts = mt_seeds.get()
    if len(mts) == 0:
        msg.showinfo("warn","At least one mod_type seed must be provided")
        return
    valid,vmsg = validate(mts,'A',1)
    if not valid:
        msg.showinfo("warn",'Mod Type: '+vmsg)
        return
    
    sts = st_seeds.get()  # Scan types are optional
    valid,vmsg = validate(sts,'A',1)
    if not valid:
        msg.showinfo("warn",'Scan Type: '+vmsg)
        return

    ############# RF ###################
    rfs = rf_seeds.get()
    sdev = rf_std.get()
    cnt = rf_sampls.get()
    
    # validate seeds
    if len(rfs) == 0:
        msg.showinfo("warn","At least one RF seed must be provided")
        return
    valid,vmsg = validate(rfs,'I',5)
    if not valid:
        msg.showinfo("warn",'RF: '+vmsg)
        return
    
    #validate sdev
    try:
        sdev = float(sdev)
    except:
        msg.showinfo("warn","Invalid RF std dev - must be real number")
        return
    try:
        cnt = int(cnt)
    except:
        msg.showinfo("warn","Invalid RF count - must be an integer")
        return
    

    seeds.append(['ELNOT',build_list(elnots,'S')])
    seeds.append(['MT',build_list(mts,'S')])
    seeds.append(['RF',build_list(rfs,'I'),[sdev,cnt]])
    if len(sts) > 0:
        seeds.append(['ST',build_list(sts,'S')])

############# PRI ###################
    pris = pri_seeds.get()
    sdev = pri_std.get()
    cnt = pri_sampls.get()
    
    # validate seeds
    if len(pris) == 0:
        msg.showinfo("warn","At least one PRI seed must be provided")
        return
    valid,vmsg = validate(pris,'I',5)
    if not valid:
        msg.showinfo("warn",'PRI: '+vmsg)
        return
    
    #validate sdev
    try:
        sdev = float(sdev)
    except:
        msg.showinfo("warn","Invalid PRI std dev - must be real number")
        return
    try:
        cnt = int(cnt)
    except:
        msg.showinfo("warn","Invalid PRI count - must be an integer")
        return

    seeds.append(['PRI',build_list(pris,'I'),[sdev,cnt]])

    ############# PD (optiona) ###################
    pds = pd_seeds.get()
    sdev = pd_std.get()
    cnt = pd_sampls.get()
    
    if len(pds)>0:
    
        # validate seeds    
        valid,vmsg = validate(pds,'I',5)
        if not valid:
            msg.showinfo("warn",'PD: '+vmsg)
            return
    
        #validate sdev & cnt
        try:
            sdev = float(sdev)
        except:
            msg.showinfo("warn","Invalid PD std dev - must be real number")
            return
        try:
            cnt = int(cnt)
        except:
            msg.showinfo("warn","Invalid PD count - must be an integer")
            return
        seeds.append(['PD',build_list(pds,'I'),[sdev,cnt]])

    ############# SP (optiona) ###################
    sps = sp_seeds.get()
    sdev = sp_std.get()
    cnt = sp_sampls.get()
    
    if len(sps)>0:
    
        # validate seeds    
        valid,vmsg = validate(sps,'I',3)
        if not valid:
            msg.showinfo("warn",'SP: '+vmsg)
            return
    
        #validate sdev & cnt
        try:
            sdev = float(sdev)
        except:
            msg.showinfo("warn","Invalid SP std dev - must be real number")
            return
        try:
            cnt = int(cnt)
        except:
            msg.showinfo("warn","Invalid SP count - must be an integer")
            return
        seeds.append(['SP',build_list(sps,'I'),[sdev,cnt]])
    
    ############# IR (optiona) ###################
    irs = ir_seeds.get()
    sdev = ir_std.get()
    cnt = ir_sampls.get()
    
    if len(irs)>0:
    
        # validate seeds    
        valid,vmsg = validate(irs,'I',3)
        if not valid:
            msg.showinfo("warn",'IR: '+vmsg)
            return
    
        #validate sdev & cnt
        try:
            sdev = float(sdev)
        except:
            msg.showinfo("warn","Invalid IR std dev - must be real number")
            return
        try:
            cnt = int(cnt)
        except:
            msg.showinfo("warn","Invalid IR count - must be an integer")
            return
        seeds.append(['IR',build_list(irs,'I'),[sdev,cnt]])
     
    ##
    # Default geo seeds (San Antonio TX & Denver, CO)
    geo_seeds = [[29.4243,-98.491],[39.744,-104.95]]
    seeds.append(['GEO',geo_seeds,[0.05]])
    print(seeds)
    int_bldr(num_ints,delete_var,True,seeds)
    refresh_elnot_list()
    msg.showinfo("info",str(num_ints)+" intercepts generated!")


def NoImp():
    msg.showinfo("tbd","Feature not yet implemented")

def dsabl():
    msg.showinfo("tbd","This feature has been disabled")

###############################################################################
#####                 WINDOWS, PANELS & WIDGETS DEFINITION                #####
###############################################################################

root = tk.Tk()
root.title("VIPR-Lite")
root.geometry("1600x1200")
leftFrame = tk.Frame(root)
rightFrame = tk.Frame(root)
topFrame = tk.Frame(root)
botFrame = tk.Frame(root)
topFrame.pack(side=tk.TOP)
botFrame.pack(side=tk.BOTTOM)
leftFrame.pack(side=tk.LEFT)
rightFrame.pack(side=tk.RIGHT)

###############################################################################
# Create/Init ELNOT Dropdown
############################
elnot_options=['ALL']

db = db_con()
cur = db.cursor()
cur. execute("select distinct elnot from intercepts order by elnot")
dat = cur.fetchall()
cur.close()
db.close()
for notation in dat:
    elnot_options.append(notation[0])

mt_options = ['ALL','D','M','S','J','B']
parm_options = ['ALL','RF','PRI','PD','SP','IR']

elnot_var = tk.StringVar(root)
mt_var = tk.StringVar(root)
parm_var = tk.StringVar(root)

if len(elnot_options)>1:
    elnot_var.set(elnot_options[1]) # default value
else:
    elnot_var.set(elnot_options[0])
mt_var.set(elnot_options[0]) # default value
parm_var.set(elnot_options[0]) # default value

elnot_pull_down = tk.OptionMenu(topFrame, elnot_var, *elnot_options)
mt_pull_down = tk.OptionMenu(topFrame, mt_var, *mt_options)
parm_pull_down = tk.OptionMenu(topFrame, parm_var, *parm_options)
elnot_pull_down.pack(side=tk.LEFT,padx=5)
mt_pull_down.pack(side=tk.LEFT,padx=5)
parm_pull_down.pack(side=tk.LEFT,padx=5)


###############################################################################
# MENU DEFINITIONS
##################

menu = tk.Menu(root)
root.config(menu=menu)

# FILE
fileMenu = tk.Menu(menu)
menu.add_cascade(label = "File",menu=fileMenu)

import_menu=tk.Menu(fileMenu)
fileMenu.add_cascade(label = "Load Intercepts...",menu=import_menu)
import_menu.add_command(label="from Database: " +srcdb, command=show_db_import)
import_menu.add_command(label="Demo (fake) data",command=show_rand_int_pnl) #display_random_ints)
import_menu.add_command(label="from a file (TBD)",command=NoImp) #refresh_elnot_list)

delete_menu=tk.Menu(fileMenu)
fileMenu.add_cascade(label = "Delete VL2 Intercepts...",menu=delete_menu)
delete_menu.add_command(label="Delete Selected ELNOT Only",command=delete_VL2_elnot_ints)
delete_menu.add_command(label="Delete ALL intercepte",command=delete_VL2_ints)
fileMenu.add_separator()

fileMenu.add_command(label = "Export VL2 PE to: " +fq_refdb+ " (disabled)",command=dsabl) #run_export)
fileMenu.add_separator()

fileMenu.add_command(label="Create/Update VL2 DB Schema",command=update_VL2_schema)
fileMenu.add_command(label="Create/Update CNF (Export) Schema: " +fq_refdb+ " (disabled)",command=dsabl) #update_refdb_schema)
fileMenu.add_separator()

fileMenu.add_command(label="Exit",command=QuitApp)

# RUN
runMenu = tk.Menu(menu)
menu.add_cascade(label="Run", menu=runMenu)
baselineMenu = tk.Menu(runMenu)
runMenu.add_cascade(label="Generate New Baseline...", menu=baselineMenu)
baselineMenu.add_command(label="Baseline Selected ELNOT",command=gen_elnot_baseline)
baselineMenu.add_command(label="Baseline All ELNOTs",command=gen_full_baseline)
runMenu.add_separator()
runMenu.add_command(label="Stage/Process Local Data (TBD)",command=None)
runMenu.add_command(label="Retrieve/Process External Data",command=show_db_proc)
runMenu.add_separator()

runMenu.add_command(label="Delete Preserved Environment",command=delete_local_env)

# VIEW
viewMenu = tk.Menu(menu)
menu.add_cascade(label="View", menu=viewMenu)
viewMenu.add_command(label="Global Clusters",command=show_clstrs) #GetGlobalClstrs)
viewMenu.add_command(label="Parameter Groups (Modes)",command=show_modes) #GetModes)
viewMenu.add_command(label="Sites",command=show_geo) #GetSites)
viewMenu.add_separator()
viewMenu.add_command(label="Intercept Histogram",command=show_int_histo_pnl)

# TUNE
tuneMenu=tk.Menu(menu)
menu.add_cascade(label="Tune", menu=tuneMenu)
tuneMenu.add_command(label="View/Update Tuning Parms",command=show_tuning_pnl)
tuneMenu.add_separator()
tuneMenu.add_command(label= "Merge Clstrs (TBD)",command=NoImp)
tuneMenu.add_command(label= "Split Clstr (TBD)",command=NoImp)
tuneMenu.add_command(label= "Delete Clstr (TBD)",command=NoImp)
tuneMenu.add_separator()
tuneMenu.add_command(label= "Reprocess Selected Clstr",command=re_clstr)
tuneMenu.add_command(label="Re-Baseline Emitter (all clstrs/modes)",command=rebaseline)

# HELP
helpMenu = tk.Menu(menu)
menu.add_cascade(label="Help",menu=helpMenu)
helpMenu.add_command(label="Documentation",command=None) #get_help)
###########
# END MENUS
###############################################################################



###############################################################################
#                      Database Query Panel
###############################################################################
dbq_pnl = tk.Frame(botFrame,bg='blue')
field_options=('','ELNOT','mod_type','country_code','collector','rd_out_stat','num_bursts')
criteria_options=('','=','<','>','!=','IN','NOT IN')

field_var1 = tk.StringVar(dbq_pnl)
field_var1.set(field_options[0])
field_var2 = tk.StringVar(dbq_pnl)
field_var2.set(field_options[0])
field_var3 = tk.StringVar(dbq_pnl)
field_var3.set(field_options[0])
field_var4 = tk.StringVar(dbq_pnl)
field_var4.set(field_options[0])
field_var5 = tk.StringVar(dbq_pnl)
field_var5.set(field_options[0])

criteria_var1 = tk.StringVar(dbq_pnl)
criteria_var1.set(criteria_options[0])

criteria_var2 = tk.StringVar(dbq_pnl)
criteria_var2.set(criteria_options[0])
criteria_var3 = tk.StringVar(dbq_pnl)
criteria_var3.set(criteria_options[0])
criteria_var4 = tk.StringVar(dbq_pnl)
criteria_var4.set(criteria_options[0])
criteria_var5 = tk.StringVar(dbq_pnl)
criteria_var5.set(criteria_options[0])

field1 = tk.OptionMenu(dbq_pnl, field_var1, *field_options)
field2 = tk.OptionMenu(dbq_pnl, field_var2, *field_options)
field3 = tk.OptionMenu(dbq_pnl, field_var3, *field_options)
field4 = tk.OptionMenu(dbq_pnl, field_var4, *field_options)
field5 = tk.OptionMenu(dbq_pnl, field_var5, *field_options)

criteria1 = tk.OptionMenu(dbq_pnl, criteria_var1, *criteria_options)
criteria2 = tk.OptionMenu(dbq_pnl, criteria_var2, *criteria_options)
criteria3 = tk.OptionMenu(dbq_pnl, criteria_var3, *criteria_options)
criteria4 = tk.OptionMenu(dbq_pnl, criteria_var4, *criteria_options)
criteria5 = tk.OptionMenu(dbq_pnl, criteria_var5, *criteria_options)


db_val1 = tk.Entry(dbq_pnl,bd=3)
db_val2 = tk.Entry(dbq_pnl,bd=3)
db_val3 = tk.Entry(dbq_pnl,bd=3)
db_val4 = tk.Entry(dbq_pnl,bd=3)
db_val5 = tk.Entry(dbq_pnl,bd=3)
min_lat = tk.Entry(dbq_pnl,bd=3)
max_lat = tk.Entry(dbq_pnl,bd=3)
min_lon = tk.Entry(dbq_pnl,bd=3)
max_lon = tk.Entry(dbq_pnl,bd=3)

db_title = tk.Label(dbq_pnl,text="Database Query Panel",bg='blue',fg='white')
db_numdays_label = tk.Label(dbq_pnl,text="Number of days*",bg='blue',fg='white')
db_numdays_entry=tk.Entry(dbq_pnl,bd=3)

db_label1 = tk.Label(dbq_pnl,text="Optional Query Criteria",bg='blue',fg='white')
db_label2 = tk.Label(dbq_pnl,text="Field",bg='blue',fg='white')
db_label3 = tk.Label(dbq_pnl,text="Criteria",bg='blue',fg='white')
db_label4 = tk.Label(dbq_pnl,text="Value",bg='blue',fg='white')
db_geo_label = tk.Label(dbq_pnl,text = "Geo Bounding Box Coordinates",bg='blue',fg='white')
min_lat_label = tk.Label(dbq_pnl,text = "Min Lat",bg='blue',fg='white')
max_lat_label = tk.Label(dbq_pnl,text = "Max Lat",bg='blue',fg='white')
min_lon_label = tk.Label(dbq_pnl,text = "Min Lon",bg='blue',fg='white')
max_lon_label = tk.Label(dbq_pnl,text = "Max Lon",bg='blue',fg='white')
db_clear_btn = tk.Button(dbq_pnl,text = "Clear All", highlightbackground='blue',fg='red',command=clear_db_criteria)
db_exit_btn = tk.Button(dbq_pnl,text='X',command=hide_db_pnls)

# import widgets
db_lagdays_entry=tk.Entry(dbq_pnl,bd=3)
db_lagdays_label=tk.Label(dbq_pnl,text = "Lag days",bg='blue',fg='white')
db_count_btn=tk.Button(dbq_pnl,text="Count",highlightbackground='blue',fg='red',command=count_vipr_ints)
db_import_btn = tk.Button(dbq_pnl,text = "Import", highlightbackground='blue',fg='red',command=import_vipr_ints)

#proc widgets
db_process_btn = tk.Button(dbq_pnl,text = "Load & Process", highlightbackground='blue',fg='red',command=proc_new_dat)

# Pack static widgets
db_exit_btn.grid(row=0,column=6,sticky='ne')
db_title.grid(row=0,column=3)
db_numdays_label.grid(row=2,column=1)
db_numdays_entry.grid(row=2,column=2)
db_label1.grid(row=5,column=3)
db_label2.grid(row=6,column=1)
db_label3.grid(row=6,column=3)
db_label4.grid(row=6,column=5)

field1.grid(row=7,column=1)
criteria1.grid(row=7,column=3)
db_val1.grid(row=7,column=5)

field2.grid(row=8,column=1)
criteria2.grid(row=8,column=3)
db_val2.grid(row=8,column=5)

field3.grid(row=9,column=1)
criteria3.grid(row=9,column=3)
db_val3.grid(row=9,column=5)

field4.grid(row=10,column=1)
criteria4.grid(row=10,column=3)
db_val4.grid(row=10,column=5)

field5.grid(row=11,column=1)
criteria5.grid(row=11,column=3)
db_val5.grid(row=11,column=5)

db_geo_label.grid(row=13,column=3)
min_lat_label.grid(row=14,column=1)
min_lat.grid(row=14,column=2)
max_lat_label.grid(row=14,column=4)
max_lat.grid(row=14,column=5)

min_lon_label.grid(row=15,column=1)
min_lon.grid(row=15,column=2)
max_lon_label.grid(row=15,column=4)
max_lon.grid(row=15,column=5)
db_clear_btn.grid(row=16,column=3)


####################
# END DB QUERY Panel
###############################################################################



###############################################################################
#                           Tuning Panel
###############################################################################
tuning_pnl = tk.Frame(botFrame,bg='blue',bd=3)
tuningLabel = tk.Label(tuning_pnl,text="Tuning Panel",bg='blue',fg='white')
tuningLabel2 = tk.Label(tuning_pnl,text=" ",bg='blue',fg='white')

update_button = tk.Button(tuning_pnl,text="Update Tuning Parms",highlightbackground='blue',fg='red',command=update_parms)
retrieve_button = tk.Button(tuning_pnl,text="Retrieve Current Parms",highlightbackground='blue',fg='red',command=get_current_parms)

#Create parm drop_down 
#pvar_options=['RF','PRI','PD','SP','IR']
#pvar = tk.StringVar(root)
#pvar.set(pvar_options[0]) # default value
#parm_pull_down = tk.OptionMenu(tuning_pnl, pvar, *pvar_options)
#w.grid(row=0,column=1)  #pack(side=tk.TOP)

h_entry = tk.Entry(tuning_pnl,bd=3)
hlabel =tk.Label(tuning_pnl,text="horizon:",bg='blue',fg='white')

t_entry = tk.Entry(tuning_pnl,bd=3)
tlabel =tk.Label(tuning_pnl,text="threshold:",bg='blue',fg='white')

p_entry = tk.Entry(tuning_pnl,bd=3)
plabel =tk.Label(tuning_pnl,text="percentile:",bg='blue',fg='white')

tuning_exit_btn=tk.Button(tuning_pnl,text='X',command=hide_db_pnls)

# Pack static content
#parm_pull_down.grid(row=2,column=1)
tuningLabel.grid(row=1,column=4)
retrieve_button.grid(row=2,column=4)
hlabel.grid(row=3,column=1)
h_entry.grid(row=3,column=2)
tlabel.grid(row=3,column=3)
t_entry.grid(row=3,column=4)
plabel.grid(row=3,column=5)
p_entry.grid(row=3,column=6)
update_button.grid(row=3,column=7)
tuningLabel2.grid(row=4,column=1)
tuning_exit_btn.grid(row=1,column=7,sticky='ne')

##################
# END Tuning Panel
###############################################################################



###############################################################################
#                  Random Intercept Generation Panel
###############################################################################
rand_int_pnl = tk.Frame(botFrame,bg='blue')
elnot_seeds = tk.Entry(rand_int_pnl,bd=3)
mt_seeds = tk.Entry(rand_int_pnl,bd=3)
st_seeds = tk.Entry(rand_int_pnl,bd=3)
rf_seeds = tk.Entry(rand_int_pnl,bd=3)
pri_seeds = tk.Entry(rand_int_pnl,bd=3)
pd_seeds = tk.Entry(rand_int_pnl,bd=3)
sp_seeds = tk.Entry(rand_int_pnl,bd=3)
ir_seeds = tk.Entry(rand_int_pnl,bd=3)
rf_std = tk.Entry(rand_int_pnl,bd=3)
pri_std = tk.Entry(rand_int_pnl,bd=3)
rf_sampls = tk.Entry(rand_int_pnl,bd=3)
pri_sampls = tk.Entry(rand_int_pnl,bd=3)
pd_std = tk.Entry(rand_int_pnl,bd=3)
pd_sampls = tk.Entry(rand_int_pnl,bd=3)
sp_std = tk.Entry(rand_int_pnl,bd=3)
sp_sampls = tk.Entry(rand_int_pnl,bd=3)
ir_std = tk.Entry(rand_int_pnl,bd=3)
ir_sampls = tk.Entry(rand_int_pnl,bd=3)

num_rand_ints = tk.Entry(rand_int_pnl,bd=3)
num_rand_label = tk.Label(rand_int_pnl,text="Number of Intercepts*",bg='blue',fg='white')

# Set some defaults
num_rand_ints.insert(0,'1000')
elnot_seeds.insert(0,'AAAAA,BBBBB')
mt_seeds.insert(0,'D,M')
rf_seeds.insert(0,'2150,6000')
rf_std.insert(0,'5')
rf_sampls.insert(0,'200')
pri_std.insert(0,'2')
pri_sampls.insert(0,'200')

pri_seeds.insert(0,'100,200,300')

    
rand_exit_btn = tk.Button(rand_int_pnl,text='X',command=hide_db_pnls)
rand_int_btn = tk.Button(rand_int_pnl,text="Generate",highlightbackground='blue',fg='red',command=gen_rand_ints) #count_vipr_ints)

rand_int_label1 = tk.Label(rand_int_pnl,text="Random Intercept Generation",bg='blue',fg='white')
rand_int_label2 = tk.Label(rand_int_pnl,text="(separate seeds by commas)",bg='blue',fg='white')
rand_std_label = tk.Label(rand_int_pnl,text="Standard Deviation",bg='blue',fg='white')
rand_sample_label = tk.Label(rand_int_pnl,text="Samples",bg='blue',fg='white')
append_var = tk.BooleanVar(value=True)
append_btn = tk.Checkbutton(rand_int_pnl, text="Append", bg='blue',fg='white',variable=append_var)

elnot_seed_label = tk.Label(rand_int_pnl,text = "ELNOT seeds*",bg='blue',fg='white')
mt_seed_label = tk.Label(rand_int_pnl,text = "MOD TYPE seeds*",bg='blue',fg='white')
st_seed_label = tk.Label(rand_int_pnl,text = "SCAN TYPE seeds",bg='blue',fg='white')
rf_seed_label = tk.Label(rand_int_pnl,text = "RF seeds*",bg='blue',fg='white')
pri_seed_label = tk.Label(rand_int_pnl,text = "PRI seeds*",bg='blue',fg='white')
pd_seed_label = tk.Label(rand_int_pnl,text = "PD seeds",bg='blue',fg='white')
sp_seed_label = tk.Label(rand_int_pnl,text = "SP seeds",bg='blue',fg='white')
ir_seed_label = tk.Label(rand_int_pnl,text = "IR seeds",bg='blue',fg='white')



# Pack static content
rand_exit_btn.grid(row=0,column=6)
rand_int_label1.grid(row=0,column=4)
rand_int_label2.grid(row=1,column=4)
rand_std_label.grid(row=5,column=4)
rand_sample_label.grid(row=5,column=5)
elnot_seed_label.grid(row=2,column=2)
elnot_seeds.grid(row=2,column=3)
mt_seed_label.grid(row=3,column=2)
mt_seeds.grid(row=3,column=3)
st_seed_label.grid(row=4,column=2)
st_seeds.grid(row=4,column=3)

rf_seed_label.grid(row=6,column=2)
rf_seeds.grid(row=6,column=3)
rf_std.grid(row=6,column=4)
rf_sampls.grid(row=6,column=5)

pri_seed_label.grid(row=7,column=2)
pri_seeds.grid(row=7,column=3)
pri_std.grid(row=7,column=4)
pri_sampls.grid(row=7,column=5)
    
pd_seed_label.grid(row=8,column=2)
pd_seeds.grid(row=8,column=3)
pd_std.grid(row=8,column=4)
pd_sampls.grid(row=8,column=5)

sp_seed_label.grid(row=9,column=2)
sp_seeds.grid(row=9,column=3)
sp_std.grid(row=9,column=4)
sp_sampls.grid(row=9,column=5)

ir_seed_label.grid(row=10,column=2)
ir_seeds.grid(row=10,column=3)
ir_std.grid(row=10,column=4)
ir_sampls.grid(row=10,column=5)

num_rand_label.grid(row=3,column=4)
append_btn.grid(row=2,column=5)
num_rand_ints.grid(row=3,column=5)
rand_int_btn.grid(row=11,column=4)
#######################################
# End Random Intercept Generation Panel
###############################################################################



###############################################################################
#                          HISTOGRAM WINDOW
###############################################################################
histo_win = tk.Frame(rightFrame)

h_subframe = tk.Frame(histo_win)
h_plot_frame=tk.Frame(h_subframe)
h_btn_frame = tk.Frame(h_subframe)

pop_histo_btn = tk.Button(h_btn_frame,text='Pop Out',command=pop_histo)
clr_histo_btn = tk.Button(h_btn_frame,text = 'Close Histo', command=clear_histo)
clr_histo_btn.pack(side=tk.LEFT)
pop_histo_btn.pack(side=tk.LEFT)
h_btn_frame.pack(side=tk.TOP)
h_plot_frame.pack()

ih_subframe = tk.Frame(histo_win)
ih_plot_frame =tk.Frame(ih_subframe)
ih_btn_frame = tk.Frame(ih_subframe)


clr_int_hist_btn = tk.Button(ih_btn_frame,text='Close Int Histo',command=clear_int_histo)
clr_int_hist_btn.pack(side=tk.LEFT)
ih_btn_frame.pack(side=tk.TOP)
ih_plot_frame.pack()


#####################
# End Histogram Panel
###############################################################################



###############################################################################
#                          Cluster Panel
###############################################################################
clstr_pnl = tk.Frame(leftFrame)
clstr_cols=('clstr_id','elnot','mod_type','parm','parm_min','parm_max')
clstr_tree = ttk.Treeview(clstr_pnl, column=clstr_cols, show='headings')
# define headings
clstr_tree.heading('clstr_id', text="CLSTR_ID")  #Note: ID by col name or #
clstr_tree.heading('elnot', text="ELNOT")
clstr_tree.heading("#3", text="PARM")
clstr_tree.heading("#4", text="MOD_TYPE")
clstr_tree.heading("#5", text="MIN")
clstr_tree.heading("#6", text="MAX")
# Center align data content
clstr_tree.column("#1", anchor=tk.CENTER,width=75)  #centers data in the column
clstr_tree.column("#2", anchor=tk.CENTER,width=75)
clstr_tree.column("#3", anchor=tk.CENTER,width=75)
clstr_tree.column("#4", anchor=tk.CENTER,width=75)
clstr_tree.column("#5", anchor=tk.CENTER,width=150)
clstr_tree.column("#6", anchor=tk.CENTER,width=150)
clstr_tree.bind('<<TreeviewSelect>>', clstr_selected)

clstr_title = tk.Label(clstr_pnl,text="Global Clusters")
clstr_title.pack(side=tk.TOP)
clstr_btn_frame=tk.Frame(clstr_pnl)
clstr_exit_btn = tk.Button(clstr_btn_frame,text='close',command=hide_clstrs)
clstr_rfrsh_btn = tk.Button(clstr_btn_frame,text='refresh',command=show_clstrs)
clstr_rfrsh_btn.pack(side = tk.LEFT)
clstr_exit_btn.pack(side=tk.LEFT)
clstr_btn_frame.pack(side=tk.BOTTOM)

###################
# END Cluster Panel
###############################################################################



##############################################################################
#                             Groups Panel
##############################################################################
grps_pnl = tk.Frame(leftFrame)
mode_cols = ('group_id','elnot','mod_type','firing_order','descr')
mode_tree = ttk.Treeview(grps_pnl, column=mode_cols, show='headings')
mode_tree.heading('group_id', text="GROUP_ID")  #Note: ID by col name or #
mode_tree.heading('elnot', text="ELNOT")
mode_tree.heading("#3", text="MOD_TYPE")
mode_tree.heading("#4", text="FIRING ORDER")
mode_tree.heading("#5", text="DESCR")
mode_tree.column("#1", anchor=tk.CENTER,width=75)  #centers data in the column
mode_tree.column("#2", anchor=tk.CENTER,width=75)
mode_tree.column("#3", anchor=tk.CENTER,width=75)
mode_tree.column("#4", anchor=tk.CENTER,width=300)
mode_tree.column("#5", width=300)

grps_title = tk.Label(grps_pnl,text="Parameter Groups (Modes)")
grps_title.pack(side=tk.TOP)
grps_btn_frame = tk.Frame(grps_pnl)
grps_rfrsh_btn = tk.Button(grps_btn_frame,text='refresh',command=show_modes)
grps_exit_btn = tk.Button(grps_btn_frame,text='close',command=hide_modes)
grps_rfrsh_btn.pack(side = tk.LEFT)
grps_exit_btn.pack(side=tk.LEFT)
grps_btn_frame.pack(side=tk.BOTTOM)
mode_tree.bind('<<TreeviewSelect>>', mode_selected)

##################
# End Groups Panel
###############################################################################


##############################################################################
#                               Parms Panel
##############################################################################
parms_pnl = tk.Frame(rightFrame)
pri_cols = ('idx','MIN','MAX')
rf_cols = ('global_id','global_min','global_max','local_min','local_max')
pri_tree  = ttk.Treeview(parms_pnl, column=pri_cols, show='headings')
rf_tree = ttk.Treeview(parms_pnl, column=rf_cols, show='headings')
pd_tree = ttk.Treeview(parms_pnl, column=rf_cols, show='headings')
sp_tree = ttk.Treeview(parms_pnl, column=rf_cols, show='headings')
ir_tree = ttk.Treeview(parms_pnl, column=rf_cols, show='headings')

pri_tree.heading('idx', text="PRI_NUM")  #Note: ID by col name or #
pri_tree.heading('MIN', text="MIN")
pri_tree.heading("#3", text="MAX")

rf_tree.heading('global_id', text="RF Clstr_ID")
rf_tree.heading('global_min', text="GLOBAL MIN")
rf_tree.heading('global_max', text="GLOBAL MAX")
rf_tree.heading('local_min', text="LOCAL MIN")
rf_tree.heading('local_max', text="LOCAL MAX")

pd_tree.heading('global_id', text="PD Clstr_ID")
pd_tree.heading('global_min', text="GLOBAL MIN")
pd_tree.heading('global_max', text="GLOBAL MAX")
pd_tree.heading('local_min', text="LOCAL MIN")
pd_tree.heading('local_max', text="LOCAL MAX")

sp_tree.heading('global_id', text="SP Clstr_ID")
sp_tree.heading('global_min', text="GLOBAL MIN")
sp_tree.heading('global_max', text="GLOBAL MAX")
sp_tree.heading('local_min', text="LOCAL MIN")
sp_tree.heading('local_max', text="LOCAL MAX")

ir_tree.heading('global_id', text="IR Clstr_ID")
ir_tree.heading('global_min', text="GLOBAL MIN")
ir_tree.heading('global_max', text="GLOBAL MAX")
ir_tree.heading('local_min', text="LOCAL MIN")
ir_tree.heading('local_max', text="LOCAL MAX")

pri_tree.column("#1", anchor=tk.CENTER,width=75)  #centers data in the column
pri_tree.column("#2", anchor=tk.CENTER,width=150)
pri_tree.column("#3", anchor=tk.CENTER,width=150)

rf_tree.column("#1", anchor=tk.CENTER,width=75)  #centers data in the column
rf_tree.column("#2", anchor=tk.CENTER,width=150)
rf_tree.column("#3", anchor=tk.CENTER,width=150)
rf_tree.column("#4", anchor=tk.CENTER,width=150)
rf_tree.column("#5", anchor=tk.CENTER,width=150)

pd_tree.column("#1", anchor=tk.CENTER,width=75)  #centers data in the column
pd_tree.column("#2", anchor=tk.CENTER,width=150)
pd_tree.column("#3", anchor=tk.CENTER,width=150)
pd_tree.column("#4", anchor=tk.CENTER,width=150)
pd_tree.column("#5", anchor=tk.CENTER,width=150)

sp_tree.column("#1", anchor=tk.CENTER,width=75)  #centers data in the column
sp_tree.column("#2", anchor=tk.CENTER,width=150)
sp_tree.column("#3", anchor=tk.CENTER,width=150)
sp_tree.column("#4", anchor=tk.CENTER,width=150)
sp_tree.column("#5", anchor=tk.CENTER,width=150)

ir_tree.column("#1", anchor=tk.CENTER,width=75)  #centers data in the column
ir_tree.column("#2", anchor=tk.CENTER,width=150)
ir_tree.column("#3", anchor=tk.CENTER,width=150)
ir_tree.column("#4", anchor=tk.CENTER,width=150)
ir_tree.column("#5", anchor=tk.CENTER,width=150)

parms_exit_btn = tk.Button(parms_pnl,text='close parms',command=hide_mode_parms)

##################
# End Parms Panel
###############################################################################



##############################################################################
#                                Geo Panel
##############################################################################
geo_pnl = tk.Frame(leftFrame)
site_cols = ('site_id','elnot','sma','smi','orient','lat','lon','heard_count','ROA')
site_tree = ttk.Treeview(geo_pnl, column=site_cols, show='headings')
site_tree.heading('site_id', text="SITE_ID")  #Note: ID by col name or #
site_tree.heading('elnot', text="ELNOT")
site_tree.heading("#3", text="SMA")
site_tree.heading("#4", text="SMI")
site_tree.heading("#5", text="ORIENT")
site_tree.heading("#6", text="LAT")
site_tree.heading("#7", text="LON")
site_tree.heading("#8", text="HEARD_COUNT")
site_tree.heading("#9", text="ROA")
site_tree.column("#1", anchor=tk.CENTER,width=75)  #centers data in the column
site_tree.column("#2", anchor=tk.CENTER,width=75)
site_tree.column("#3", anchor=tk.CENTER,width=100)
site_tree.column("#4", anchor=tk.CENTER,width=100)
site_tree.column("#5", anchor=tk.CENTER,width=75)
site_tree.column("#6", anchor=tk.CENTER,width=125)
site_tree.column("#7", anchor=tk.CENTER,width=125)
site_tree.column("#8", anchor=tk.CENTER,width=100)
site_tree.column("#9", anchor=tk.CENTER,width=75)
#site_tree.bind('<<TreeviewSelect>>', site_selected)

geo_title = tk.Label(geo_pnl,text="Geolocations (Sites)")
geo_title.pack(side=tk.TOP)
geo_btn_frame = tk.Frame(geo_pnl)
geo_rfrsh_btn=tk.Button(geo_btn_frame,text='refresh',command=show_geo)
geo_exit_btn=tk.Button(geo_btn_frame,text='close geo',command=hide_geo)
geo_rfrsh_btn.pack(side=tk.LEFT)
geo_exit_btn.pack(side=tk.LEFT)
geo_btn_frame.pack(side=tk.BOTTOM)

###############
# END Geo Panel
###############################################################################

##############################################################################
#                         int_histo_panel
##############################################################################
int_hist_pnl = tk.Frame(botFrame,bg='blue')
ih_title = tk.Label(int_hist_pnl,text="Use pull-downs to select parms",bg='blue',fg='white')
ih_min_label = tk.Label(int_hist_pnl,text="Min:",bg='blue',fg='white')
ih_max_label = tk.Label(int_hist_pnl,text="Max:",bg='blue',fg='white')
ih_min_entry = tk.Entry(int_hist_pnl,bd=3)
ih_max_entry = tk.Entry(int_hist_pnl,bd=3)
ih_plot_btn = tk.Button(int_hist_pnl,text = "Plot", highlightbackground='blue',fg='red',command=plot_int_histo)
ih_exit_btn = tk.Button(int_hist_pnl,text='X',command = hide_db_pnls)

# Pack static content
ih_title.grid(row=0,column=2)
ih_min_label.grid(row=1,column=0)
ih_min_entry.grid(row=1,column=1)
ih_max_label.grid(row=1,column=3)
ih_max_entry.grid(row=1,column=4)
ih_plot_btn.grid(row=2,column=2)
ih_exit_btn.grid(row=0,column=5)



    

show_clstrs()
show_modes()
show_geo()
root.mainloop()
