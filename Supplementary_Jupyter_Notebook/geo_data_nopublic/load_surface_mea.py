#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 10:11:51 2020

@author: brennanbrunsvik
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fname = 'surface_meas_boncio'

sm = pd.read_table(fname)
no_meas = np.isnan(sm.Lat)
sm = sm[~no_meas]
nostrike = sm.strike == '-'
sm.strike[nostrike] = np.nan
sm.strike = sm.strike.astype(float)

# plt.scatter(sm.Long, sm.Lat)
codes = sm.code
pag = np.zeros(codes.shape, dtype = bool)

codes = codes.to_numpy()
for icode, code in enumerate(codes):
    print(code)
    if code[0] == 'P':
        pag[icode] = True
        
sm = sm.strike.to_numpy()
pagstrike = sm[pag]

# plt.hist(sm.strike)
# plt.show()

plt.scatter(sm.Long, sm.Lat, c = sm.strike, vmin = 90, vmax = 160)
plt.colorbar()