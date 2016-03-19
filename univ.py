
r1 = 2.2; r2 = 6
corr = .00001

import math
import MDAnalysis
import pylab
import numpy as np
from numpy.linalg import norm
import matplotlib as mlab
from datetime import datetime, timedelta

def normalV(A, B, C, direct):
    vec = np.array([0, 0, 0])
    u = A - B
    v = B - C
    vec[0] = u[1]*v[2]-u[2]*v[1]
    vec[1] = u[2]*v[0]-u[0]*v[2]
    vec[2] = u[0]*v[1]-u[1]*v[0]
    #vec = vec/norm(vec)
    if(np.dot(vec, direct) < 0):
        vec *= -1
    return vec

def angleOH(xO, xH, nSurf):
    d   = xH - xO; #h2 = h2 - pos
    cos = np.dot(nSurf, d)/norm(nSurf)/norm(d)
    #cos2 = np.dot(n, h2)/norm(n)/norm(h2)
    #cos  = cos1
    #if cos2 < cos1: cos = cos2 
    return cos
    
u = MDAnalysis.Universe("mdconf.gro","trj.trr")

g_part = u.select_atoms("name Ti or name O")

N_part = len(g_part)
g_sol  = u.residues[N_part:]
g_o    = g_sol.select_atoms("name OW")
N_sol  = len(g_sol)
print "N_sol = %d, N_part = %d" % (N_sol, N_part)

str_arg = "around " + str(r1) + " (name Ti or name O)"
first  = u.select_atoms(str_arg)
first  = first.select_atoms("name OW")
#print first.names
print "first N = %d" % (len(first))

xlist = []; zlist = []#; num = 0
now = datetime.now()
print(now)
for ts in u.trajectory:
    str_arg = "around " + str(r1) + " (name Ti or name O)"
    first  = u.select_atoms(str_arg).select_atoms("name OW")
    #first  = first.select_atoms("name OW")
    print ts.frame
    if(ts.frame == 702):
        print ts.frame
        break
    for i in range(0, len(first)):
        ind = first.indices[i] + 1
        str_arg = "around 4 (bynum " + str(ind) + " )"
        sel = u.select_atoms(str_arg)
        sel = sel.select_atoms("name Ti or name O")
            
        pos  = first.positions[i]
        near = sel.positions
        d = near - pos
        d = norm(d, axis=1)
        d = np.sort(d)
        d = d[0:3]
        str_arg = "around " + str(np.max(d)+corr) + " (bynum " + str(ind) + " )"
        sel = u.select_atoms(str_arg).select_atoms("name Ti or name O")
        
        a = sel.positions[0]
        b = sel.positions[1]
        c = sel.positions[2]
        n = normalV(a, b, c, pos)
        
        h1 = u.atoms.positions[ind]        
        zlist.append(angleOH(pos, h1, n))
    #num += 1
    #print "%(i)d: %(cos)f %(sel)r" % vars()
now = datetime.now()
print(now)    
ndata = len(zlist); zlist = np.sort(zlist)
ylist = np.zeros(1000); num = 0
for i in range(0, 1000):
    dmax = -1 + (i+1)*.002
    while zlist[num] < dmax:
        ylist[i] += 1
        num += 1
        if num == ndata: break
    if num == ndata: break

xlist = frange(-1,.998,0.002)
pylab.plot (xlist, ylist, "gD")
pylab.show()
#for ts in u.trajectory:
#    for i in range(0,10):
#        pnt = g_o[i]
#        near  = u.select_atoms("around 0.3 pnt")        
#        #print len(near)

#z = []; x = []; dz = [] 
#for r in frange(2., 7., .02):
#    str_arg = "around " + str(r) + " (name Ti or name O)"
#    first  = u.select_atoms(str_arg)
#    first  = first.select_atoms("name OW")
#    z.append(len(first))
#    x.append(r)
#
#for i in range(0,len(z)):
#    val = z[i]*()/(2+i*.02)**2
#    dz.append(30*val)

   