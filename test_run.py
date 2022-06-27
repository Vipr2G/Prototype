#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 08:42:24 2021

@author: russelmiller
"""

from VL2bb import build_baseline
from VL2residue_miner import new_residue_miner
from VL2associator import new_associator
from VL2geo_associator import associate_geo
from VL2geominer import proc_geo
from VL2db_utils import delete_env, pe_stats
from VL2test_utils import int_bldr
from VL2pretty import pretty_print

print()
print("##### Test Run: Deleting IPE")
print()
delete_env(0)

print()
print("##### Test Run: Creating baseline intercepts")
print()
int_bldr(2500)

print()
print("##### Test Run: Building Baseline")
print()
build_baseline('BBBBB')

print()
print("##### Test Run: Creating new intercepts for association")
print()
int_bldr(2500,False)

print()
print("##### Test Run: Running parametric association")
print()
new_associator()

print()
print("##### Test Run: Running geo association")
print()
associate_geo()

print()
print("##### Test Run: Mining parametric residue")
print()
new_residue_miner()

print()
print("##### Test Run: Mining geo residue")
print()
proc_geo('BBBBB')

print()
print("##### Test Run: Printing Results")
print()
pretty_print()
pe_stats()


