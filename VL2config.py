#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 15:06:03 2021

@author: russelmiller
"""
#######################################################################################
# GUI Configuration
#######################################################################################
srcdb = 'VIPR'     # Schema qualifier for database used to load intercepts into local DB
fq_refdb = 'ref_modes'        # Schema qualifier to build/export to CNF schema (consistent with ref_db credentials)


#######################################################################################
# database Query throttling
#######################################################################################
max_db_trans = 500          # max number of query executions between commits
max_db_list_len = 100       # max length of a list for an "in" clause

#######################################################################################
# T & E Options
#######################################################################################
fq_idb = ''   # string to fully qualify selects from intercepts/intercept_pris for VL2 processing, e.g., 'vipr.' or '' for local schema
save_children = False   # True --> explicitly store child sequences, False --> just associate "child intercepts" to appropriate parent
#multimode_assoc = True  # allows intercepts to associate to all matched modes   # not currently used
local_rf = True         # include local RF clstrs on modes for compare/contrast with the global RF clusters
delete_pri_orphans = True
delete_rf_orphans = False


#######################################################################################
#  Parametric Extracrion Thresholds
#######################################################################################
min_residue_ints = 5   # minimum number of residue intercepts required for residue_miner processing
min_seq_cnt = 5        # minimum number of occurrences to extract a sequence from residue
min_mode_parm_cnt = 1  # minimum observed occurrences to attach a parm_clstr or scan type to a mode (residue miner)
min_baseline_ints = 10 # minimum number of available intercepts to run baseline_builder
min_clstr_ints = 10    # minimum number of intercepts to search for more clusters (baseline_builder)
min_clstr_peak = 4     # minimum number of parm values to define a meaningful peak during clustering
min_clstr_volume = 5   # minimum number of values in a new clstr required to preserve it
min_clstr_density = 1  # minimum density (count/range) ... set to 0 for no density constraint
min_st_count = 5       # minimum number of occurences to add a new scan_type during association


#######################################################################################
#  Clstr Housekeeping
#######################################################################################
clstr_split_thresh = 0.05   # proportion of max count to trigger a split


########################################################################################
# GeoMiner Configuration Parameters
########################################################################################
baseline_geo = True     # when True, baseline_builder will include geo_processing
gb_iters = 10           # max number of geo baseline iterations
def_roa = 5            # default radius of association assigned new sites on db insert: None, or desired ROA (nmi)
max_SMA = 5            # maximum SMA in nmi that will be processed
def_grid_res = 9.26    # grid resolution for heat man generation  (9.26km = 5nmi)
axis_min_nm = 0.5      # minimum SMA/SMI...need to investigate
min_temp_pct = 0.075   # percentage (decimal) of maximum temp to recognize a hot spot in heat hunter; Default = 0.1
min_geo_corr_ints = 5  # minimum intercept count to recognize a hot spot in heat hunter (max of pct/int count)


########################################################################################
# Geo Association Configuration Parameters
########################################################################################
LOC = 0.997          # Level of confidence expressed as a decimal
track_deltas = False # True --> Store the evolution of the sites in site_history for plotting
save_tracks = False  # True --> traditional correlator mode; False --> VIPR mode: send unassociated ints to residue
w = 0.01             #KF noise parm
hc_weighting = True  # True --> weights site sma/smi by inverse of heard_count for KF update
hc_max = 100         #   - Maximum heard_count used to weight site location
use_sqrt = False     #   - Use sqrt(heard_count) as a weighting parameter
stage_geo = True     # When True, parametric associator will automatically stage data for geo association
                     #   - otherwise, data must be manually staged in comparison_geo_snapshot prior to running geo_associator


########################################################################################
#   Default Clustering parameters (used when no entry exists in table direct_tune_parms)
########################################################################################

def_incr = 0.1  # for all parameters

#Radio Frequency (RF):
def_rf_horizon = 5.0   # (MHz)
def_rf_thresh  = 0.5
def_rf_p       = 0.95

# Pulse Repetition Interval (PRI):
def_pri_horizon = 2.0   # (usec)
def_pri_thresh  = 0.2
def_pri_p       = 0.95

# Pulse Duration (PD)
def_pd_horizon = 1.0   # (usec)
def_pd_thresh  = 0.2
def_pd_p       = 0.95

# Scan Period (SP):
def_sp_horizon = 1.0   # (sec)
def_sp_thresh  = 0.2
def_sp_p       = 0.95

# Illumination Rate (IR):
def_ir_horizon = 1.0   # (usec)
def_ir_thresh  = 0.2
def_ir_p       = 0.95


