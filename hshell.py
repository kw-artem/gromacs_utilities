import math
import MDAnalysis
import pylab
import numpy as np
from numpy.linalg import norm
import matplotlib as mlab
#import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from matplotlib.ticker import ScalarFormatter
from nmsurface import sliceSurfaceAtoms
from scipy.optimize import curve_fit
import biexponentialf as biexp

start = 0
step  = 5
end   = 25001
dtau  = .02
fqExc = 500
fqVw  = 2000
mpth = "/home/tema/mdMain/"; who = "RUT/35/"
path = mpth + who
top = path+"md/rs.gro"
trj = path + "/md/_traj.trr"
#top0 = "/home/tema/mdMain/RUT/25/md/_rs.gro"
#trj0 = "/home/tema/mdMain/RUT/25/md/_traj.trr"
r1 = 0; r2 = 2.3; D = 35.5

def surfaceGroups(indices):
    amountTi = len(np.argwhere(np.array(indices)<N_part/3))
    amountO  = N_part - amountTi
    return amountTi, amountO

rsExc = start % fqExc
rsVw  = start % fqVw
nframes = int((end - start)/step) + 1
u = MDAnalysis.Universe(top, trj)
#u = MDAnalysis.Universe(top,trj)

particle = u.select_atoms("name Ti or name O")
woxygens = u.select_atoms("resname SOL and name OW")

N_part = len(particle)
N_oxy  = len(woxygens)
print "N_oxy = %d, N_part = %d" % (N_oxy, N_part)
shellstate = np.zeros((nframes, N_oxy), np.bool)
fnmolc = open(path+"calc/rsT_nmolc.txt","w")
fnmolc.write('@title "Number of water molecules vs $\delta$ time "\n')
fnmolc.write('@subtitle "within a layer {0} <= r < {1} A"\n'.format(r1, r2))
nmolc = np.zeros(nframes, np.int32)
fsgrs   = open(path+"calc/numOfGrsOnSurf.txt","a")
print(datetime.now())
k = 0
for ts in u.trajectory[start:end:step]:
    frame = ts.frame
    if(frame % fqVw == rsVw):    
        print ts.frame
    if(frame % fqExc == rsExc):
        coordPart = particle.positions
        indices   = sliceSurfaceAtoms(coordPart, D-3)
        qsTi, qsO = surfaceGroups(indices)
        fsgrs.write("{0} {1} ".format(qsTi, qsO))
    coordSurf = particle.positions[indices]
    coordOxy  = woxygens.positions
    currstate = np.zeros(N_oxy, dtype = np.int32)
    i = 0
    for pos in coordOxy:
        dx = coordSurf - pos
        dl = norm(dx, axis=1)
        arg  = np.argwhere(dl<r2)
        arg  = np.delete(arg, np.argwhere(dl[arg]<r1))
        arg  = arg.reshape((arg.shape[0],))
        narg = len(arg)
        nmolc[k] += narg
        if narg > 0: shellstate[k, i] = 1
        i += 1
    fnmolc.write("{0:14.6f}{1:14.6f}\n".format(k*dtau*step, nmolc[k]))
    k += 1
fnmolc.close()
fsgrs.close()
print(datetime.now())
#%%
fout = open(path+"calc/shellstate.txt","w")
for k in range(nframes):
    for i in range(N_oxy):
        fout.write("{} ".format(int(shellstate[k, i])))
    fout.write("\n")
fout.close()
'''
#%%
total = 3000; fromt = 0
corrfunc = np.zeros(total)
rtcf  = np.zeros(total, dtype = np.float64)
for init in range(nframes)[fromt:(1-total)]:
    k = 0
    for i in range(nframes)[init:(init+total)]:
        corrfunc[k] = sum(shellstate[init]*shellstate[i])
        k += 1
    rtcf += corrfunc/nmolc[init]
rtcf /= (nframes - total - fromt + 1)
#%%
fcfunc = open("rsT_corFunc.txt","w")
for i in range(total):
    fcfunc.write("{0:10.6f}{1:10.6f}\n".format(i*dtau*step, rtcf[i]))
fnmolc.close()
#%%
print 'done'
#%%
def func(x, A, atau, btau):
    return A*np.exp(-x/atau)+(1-A)*np.exp(-x/btau)
xdata = np.arange(total)[10:1000]
ydata = rtcf[10:1000]
popt, pcov = curve_fit(func,xdata,ydata)
perr = np.sqrt(np.diag(pcov))

#a = np.array([.5, .5, 5000, 150000])
#%%
k=0; err=[]
#while k<30:#not stConverg:
#    a_prev = np.append(a_prev, a)
#    grad = biexp.gradient(x, y, a)
#    mtx2 = biexp.secondDerivativeMatrix(x, y, a)
#    s = -grad/norm(grad)
#    l = -grad.dot(s)/(s.dot(mtx2.dot(s)))
#    a+= l*s
#    err.append(abs(a-a_prev[k])/a_prev[k])
#    k+=1
'''
