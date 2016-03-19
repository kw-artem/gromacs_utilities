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

step  = 5
dtau  = .02
r1 = 0; r2 = 2.3
mpth = "/home/tema/mdMain/"; who = "ANT/33/"
path = mpth + who

def evaluateCF(total, fromt=0):
    corrfunc = np.zeros(total)
    rtcf  = np.zeros(total, dtype = np.float64)
    for init in range(nframes)[fromt:(1-total)]:
        k = 0
        for i in range(nframes)[init:(init+total)]:
            corrfunc[k] = sum(states[init]*states[i])
            k += 1
        rtcf += corrfunc/sum(states[init])
    rtcf /= (nframes - total - fromt + 1)
    return rtcf
def save(rtcf):
    fcfunc = open(path+"calc/RTCF.txt","w")
    fcfunc.write('@title "The residence time correlation function"\n')
    fcfunc.write('@subtitle "within a layer {0} <= r < {1} A"\n'.format(r1, r2))
    for i in range(len(rtcf)):
        fcfunc.write("{0:14.6f}{1:14.6f}\n".format(i*dtau*step, rtcf[i]))
    fcfunc.close()
    return None
def func(x, A, atau, btau):
    return A*np.exp(-x/atau)+(1-A)*np.exp(-x/btau)

#%%
finp = open(path+"calc/shellstate.txt","r")
s = finp.readline()
s = s.split(" ")[:-1]
N = len(s)
finp.seek(0)
nframes = 0
for s in finp.xreadlines():
    nframes += 1
states = np.zeros((nframes, N), np.bool)
finp.seek(0)
k = 0
for inp in finp.xreadlines():
    state     = np.array(inp.split(" ")[:-1], dtype = np.int16)
    states[k] = state == 1
    k += 1
finp.close()
if sum(states[-1]==1)==0:
    states=np.delete(states, nframes-1,axis=0)
    nframes = len(states)
print 'done!'
#%%
#ydata = rf275s
#xdata = np.arange(len(ydata))
#popt, pcov = curve_fit(func,xdata,ydata)
#perr = np.sqrt(np.diag(pcov))
#perr_= max(abs(func(xdata, popt[0], popt[1], popt[2])-ydata))
#%%

