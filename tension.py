# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 15:14:38 2016

@author: tema
"""

import math
import MDAnalysis
import pylab
import numpy as np
import matplotlib as mlab
from numpy.linalg import norm
from datetime import datetime, timedelta

u  = MDAnalysis.Universe("md0.gro","trj.trr")
m1 = .47867; m2 = .159994
particle = u.select_atoms("name Ti or name O")
N1 = N/3; N2 = N - N1
M1 = m1*np.ones(N1, dtype=np.float32)
M2 = m2*np.ones(N2, dtype=np.float32)
M  = np.concatenate((M1, M2), axis=0)
N = particle.n_atoms; V = 20; k = 0
ylist = []
for ts in u.trajectory:
    if(ts.frame == 2000): break
    if(ts.frame >= 1000):
        if(ts.frame % 100 == 0):    
            print ts.frame
        k += 1
        x = particle.positions
        v = particle.velocities
        f = particle.forces
        P = np.zeros((3,3), dtype = np.float32)
        
        for i in (0, 1, 2):
            for j in (0, 1, 2):
                P[i,j] = P[i,j] + ((sum(m*v[:,i]*v[:,j]))+sum(x[:,i]*f[:,j]))/V
        ylist.append((P/k).trace())

xlist = range(k); l = len(xlist)
xlist = xlist[l/2:]; ylist = ylist[l/2:];
pylab.plot (xlist, ylist, "gD")
pylab.show()
print '   done!'
