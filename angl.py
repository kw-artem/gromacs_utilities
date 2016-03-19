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

start = 0
step  = 10
end   = 25001
fqExc = 1000
fqAcc = 500
fqVw  = 2000

#top = "/home/tema/mdMain/ANT/22/md/_as.gro"
#trj = "/home/tema/mdMain/ANT/22/md/_traj.trr"
#top0 = "/home/tema/mdMain/RUT/19/md/_rs.gro"
#trj0 = "/home/tema/mdMain/RUT/19/md/_traj.trr"
mpth = "/home/tema/mdMain/"; who = "ANT/33/"
path = mpth + who
top = path+"md/as.gro"
trj = path + "/md/_traj.trr"
chkFound = False
r1 = 0; r2 = 2.3; D = 33.1

rsExc = start % fqExc
rsAcc = start % fqAcc
rsVw  = start % fqVw
def plottingHere():
    return None
def updatePickedOutWaterOxygens():
    global nsel
    global selPoints
    global listOxy
    selPoints = []; listOxy = []
    nsel = np.array([], dtype = np.int32)
    dev1 = 1; dev2 = 1.5
    coordSurf = particle.positions[indices]
    coordSol  = solution.positions
    coordOxy  = coordSol[::3]
    i = 0
    for point in coordOxy:
        dl   = norm(coordSurf - point, axis=1)
        arg  = np.argwhere(dl<r2+dev2)
        arg  = np.delete(arg, np.argwhere(dl[arg]<r1-dev1))
        arg  = arg.reshape((arg.shape[0],))
        narg = len(arg)
        if narg > 3:
            nsel = np.append(nsel, narg)
            selPoints.append(arg.tolist()) 
            listOxy.append(i)
        i += 1
    return None
def toHistogram(array, ax, bx, Np, fdest):
    ylist = np.zeros(Np); k = 0
    ndata = len(array)
    elist = np.sort(array)
    dx = (bx - ax)/Np
    fout = open(fdest, "w")
    fout.write("#{0:10d}\n".format(ndata))
    lside = ax + dx/2; rside = bx - dx/2
    for i in range(0, Np):
        dmax = ax + (i+1)*dx
        while k < ndata and elist[k] < dmax:
            ylist[i] += 1
            k += 1
        ylist[i] /= ndata
        fout.write("{0:10.6f}{1:10.6f}\n".format(dmax-dx/2, ylist[i]))
    fout.close()
    xlist = frange(lside,rside,dx)
    return xlist, ylist
def normalV(A, B, C=0., direct=None):
    vec = np.array([0., 0., 0.])
    u = A - C
    v = B - C
    vec[0] = u[1]*v[2]-u[2]*v[1]
    vec[1] = u[2]*v[0]-u[0]*v[2]
    vec[2] = u[0]*v[1]-u[1]*v[0]
    if direct != None and np.dot(vec, direct) < 0:
        vec *= -1
    return vec
def angleOH(xOH1, nSurf):
    cos = np.dot(nSurf, xOH1)/norm(nSurf)/norm(xOH1)
    return cos
def angleHH(xHH, nSurf):
    cos = np.dot(nSurf, xHH)/norm(nSurf)/norm(xHH)
    return abs(cos)
def angleDp(xOH1, xHH, nSurf):
    dp  = xOH1 + xHH/2
    cos = np.dot(nSurf, dp)/norm(nSurf)/norm(dp)
    return cos
def angleNW(xOH1, xHH, nSurf):
    nw  = normalV(xOH1, np.array([1,0,0]))
    cos = np.dot(nSurf, nw)/norm(nSurf)/norm(nw)
    return abs(cos)
def seekRound(selected):
    d = coordSurf[selected] - pos
    r = norm(d, axis=1)
    global chkFound
    chkFound = True
    tmpind = np.argsort(r)
    tri    = tmpind[:3]
    tmpdst = r[tri]
    l = tmpdst[0]
    if l < r1 or l >= r2:
            chkFound = False
#    for l in tmpdst:
#        if l < r1 or l >= r2:
#            chkFound = False
#            break
    return coordSurf[np.array(selected)[tri]]
def surfaceGroups(indices):
    amountTi = len(np.argwhere(np.array(indices)<N_part/3))
    amountO  = N_part - amountTi
    return amountTi, amountO
   
u = MDAnalysis.Universe(top, trj)

particle = u.select_atoms("name Ti or name O")
solution = u.select_atoms("resname SOL")

N_part = len(particle)
N_sol  = len(solution)
N_oxy  = N_sol/3
print "N_sol = %d, N_part = %d" % (N_sol, N_part)
print(datetime.now())

costheta = np.array([], dtype = np.float32)
cosphi   = np.array([], dtype = np.float32)
#fsgrs   = open(path+"calc/numOfGrsOnSurf.txt","a")
fout = open(path+"calc/bivardistr1.txt","w")
for ts in u.trajectory[start:end:step]:
    frame = ts.frame
    if(frame % fqVw == rsVw):    
        print ts.frame
    if(frame % fqExc == rsExc):
        coordPart = particle.positions
        indices   = sliceSurfaceAtoms(coordPart, D-3)
        nsurf = len(indices)
        #qsTi, qsO = surfaceGroups(indices)
        #fsgrs.write("{0} {1} ".format(qsTi, qsO))
        updatePickedOutWaterOxygens()
    coordSurf = particle.positions[indices]
    coordSol  = solution.positions
    coordOxy  = coordSol[::3]
    i = 0
    for ioxy in listOxy:
        pos     = coordSol[ioxy*3]
        a, b, c = seekRound(selPoints[i])
        i += 1
        if chkFound:
            h1  = coordSol[ioxy*3+1]
            h2  = coordSol[ioxy*3+2]
            gc  = (a + b + c)/3
            n   = normalV(a, b, c, pos-gc)
            oh1 = h1 - pos; hh = h2 - h1
            dip = oh1 + hh/2
            nw  = normalV(hh, dip)
            tht = np.dot(dip, n)/norm(n)/norm(dip)
            phi = np.dot(n, nw)/norm(nw)/norm(n)/sqrt(1-tht*tht)
            if phi < 0: phi *= -1
            costheta = np.append(costheta, tht)
            cosphi   = np.append(cosphi, phi)
    if(frame % fqAcc == rsAcc) or (frame == end):
        for i in xrange(len(cosphi)):
            fout.write("{0:12.8f}{1:12.8f}\n".format(cosphi[i], costheta[i]))
        costheta = np.array([], dtype = np.float32)
        cosphi   = np.array([], dtype = np.float32)
fout.close()            
print(datetime.now())

'''
# Here, one might output obtained above distributions 
  in a single .svg picture and four .txt chart-files:

#degreeFormatter = matplotlib.ticker.FormatStrFormatter(r"%g$^\circ$")
fig = plt.figure(figsize=(8,8))

x, y = toHistogram(dplist, -1., 1., 1000, path+"calc/dplist3.txt")
ax1 = fig.add_subplot(221)
ax1.plot(x, 100*y, '*', lw=3)
#ax1.plot(time, LID, 'r-', lw=2, label=r"$\theta_{\mathrm{LID}}$")
ax1.set_xlabel(r"cos $\alpha$")
ax1.set_ylabel(r"Probability, %")
#ax1.yaxis.set_major_formatter(degreeFormatter)
ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))


x, y = toHistogram(nwlist, 0., 1., 1000, path+"calc/nwlist3.txt")
ax2 = fig.add_subplot(222)
ax2.plot(x, 100*y, 'o', lw=3)
ax2.set_xlabel(r"cos $\beta$")
ax2.set_ylabel(r"Probability, %")
ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")

x, y = toHistogram(HHlist, 0., 1., 1000, path+"calc/hhlist3.txt")
ax3 = fig.add_subplot(223)
ax3.plot(x, 100*y, 'x', lw=3)
ax3.set_xlabel(r"cos $\gamma$")
ax3.set_ylabel(r"Probability, %")
ax3.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))

x, y = toHistogram(OHlist, -1., 1., 1000, path+"calc/ohlist3.txt")
ax4 = fig.add_subplot(224)
ax4.plot(x, 100*y, 'x', lw=3)
ax4.set_xlabel(r"cos $\chi$")
ax4.set_ylabel(r"Probability, %")
ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
ax4.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))


fig.subplots_adjust(left=0.12, right=0.88, bottom=0.2, wspace=0.15)

for ext in ('svg',):
    fig.savefig(path+"calc/anglesdistr.{0}".format(ext))
'''