#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 10:02:48 2020

@author: brennanbrunsvik
"""

import numpy as np
from scipy import interpolate

N = 3; M = 3
a = np.arange(N*M).reshape((N, M)) 
b = a.copy()
c = a.copy()
spl = interpolate.SmoothBivariateSpline(a, b, c)

#spl.degrees