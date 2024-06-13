#!/usr/bin/env python3
import os
import sys
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'dependencies','pyWopwop'))
import wopwop
#%%

cases_directory = os.path.dirname(__file__)
cases = 'cases.nam'
f1 = lambda a: wopwop.extract_wopwop_quant(case_directory=a, prefix = 'pressure')
f2 = lambda a: wopwop.extract_wopwop_quant(case_directory=a, prefix = 'spl_spectrum')

wopwop.apply_to_namelist([f1], cases_directory=cases_directory, cases=cases)


