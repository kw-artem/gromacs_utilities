
import math
import MDAnalysis
import pylab
import numpy as np
from numpy.linalg import norm
import matplotlib as mlab
from datetime import datetime, timedelta
from nmsurface import sliceSurfaceAtoms

start = 0
step  = 10
end   = 25001
fqExc = 500
fqAcc = 500
fqVw  = 400
nbins = 1000
ax    = 0
bx    = 6.
cutx  = 6.
dx = (bx - ax)/nbins
dev = 2; partSize = 35.5

rdf_O = np.zeros(nbins, dtype = np.float64)
rdf_H = np.zeros(nbins, dtype = np.float64)

rsExc = start % fqExc
rsAcc = start % fqAcc
rsVw  = start % fqVw
def updateWaterExceptions():
    global cons_WO
    global cons_WH
    coord_WO  = water_O.positions
    coord_WH  = water_H.positions
    coordPart = particle.positions
    averPoint = coordPart.mean(axis=0)
    dr_WO = coord_WO - averPoint
    dr_WH = coord_WH - averPoint
    dl_WO = norm(dr_WO, axis=1)
    dl_WH = norm(dr_WH, axis=1)
    cons_WO = np.argwhere(dl_WO < partSize/2 + cutx + dev)
    cons_WH = np.argwhere(dl_WH < partSize/2 + cutx + dev)
    cons_WO = cons_WO.reshape((cons_WO.shape[0],))
    cons_WH = cons_WH.reshape((cons_WH.shape[0],))
    return None
def accumulateToHistogram(portion, heap):
    ndata = len(portion)
    k = 0
    for ibin in range(nbins):
        xmax = ax + (ibin+1)*dx
        while k < ndata and portion[k] < xmax:
            heap[ibin] += 1
            k += 1
    return None
def normalizeHistogram(heap, npairs):
    divisor = 4*pi*(npairs/V)*dx
    for ibin in range(nbins):
        x = ax + (ibin+1)*dx
        heap[ibin] /= x*x*divisor
    heap /= int((end - start)/step)
    return None
rdf_O = np.zeros(nbins, dtype = np.float64)
rdf_H = np.zeros(nbins, dtype = np.float64)
def surfaceGroups(indices):
    amountTi = len(np.argwhere(np.array(indices)<N_part/3))
    amountO  = N_part - amountTi
    return amountTi, amountO

u = MDAnalysis.Universe("/home/tema/mdMain/RUT/35/md/rs.gro",
                        "/home/tema/mdMain/RUT/35/md/_traj.trr")
if end < 0: end += u.trajectory.n_frames + 1
#fids = open("_idsSurf.txt", "r")
#str  = fids.read()
#indices = [int(element) for element in str.split(' ')[1:]]

particle = u.select_atoms("name Ti or name O")
solution = u.select_atoms("resname SOL")
water_O  = solution.select_atoms("name OW")
water_H  = solution.select_atoms("name HW*")

N_part = len(particle)
N_sol  = len(solution)
N_oxy  = N_sol/3; N_hyd = N_sol - N_oxy
print "N_sol = %d, N_part = %d" % (N_sol, N_part)

distr  = []
curVolume = np.array([])
cons_WO = np.array([], dtype = np.int32)
cons_WH = np.array([], dtype = np.int32)
dl_WO   = np.array([], dtype = np.float32)
dl_WH   = np.array([], dtype = np.float32)
fsurf   = open("/home/tema/mdMain/RUT/35/calc/numOfAtsOnSurf.txt","a")
fsgrs   = open("/home/tema/mdMain/RUT/35/calc/numOfGrsOnSurf.txt","a")
print(datetime.now())
for ts in u.trajectory[start:end:step]:
    frame = ts.frame
    #if(frame == end + 1): break
    if(frame % fqVw == rsVw):    
        print ts.frame
    if(frame % fqExc == rsExc):
        coordPart = particle.positions
        indices   = sliceSurfaceAtoms(coordPart, 32.5)
        qsTi, qsO = surfaceGroups(indices)
        nsurf = len(indices)
        fsurf.write(" {}".format(nsurf))
        fsgrs.write("{0} {1} ".format(qsTi, qsO))
        updateWaterExceptions()
    coord_WO  = water_O.positions[cons_WO]
    coord_WH  = water_H.positions[cons_WH]
    coordSurf = particle.positions[indices]
    curVolume = np.append(curVolume, ts.volume)
    for xSurf in coordSurf:
        dr_WO = coord_WO - xSurf
        dr_WH = coord_WH - xSurf
        tmpWO = norm(dr_WO, axis=1)
        tmpWH = norm(dr_WH, axis=1)
        tmpWO = np.delete(tmpWO, np.argwhere(tmpWO > cutx))
        tmpWH = np.delete(tmpWH, np.argwhere(tmpWH > cutx))
        dl_WO = np.concatenate((dl_WO, tmpWO))
        dl_WH = np.concatenate((dl_WH, tmpWH))
    if(frame % fqAcc == rsAcc) or (frame == end):
        dl_WO.sort(); dl_WH.sort()
        accumulateToHistogram(dl_WO, rdf_O)
        accumulateToHistogram(dl_WH, rdf_H)
        dl_WO = np.array([], dtype = np.float32)
        dl_WH = np.array([], dtype = np.float32)
print(datetime.now())

fsurf.close()
fsgrs.close()
V = curVolume.mean()
normalizeHistogram(rdf_O, nsurf*N_oxy)
normalizeHistogram(rdf_H, nsurf*N_hyd)

f_rdf_O = open("/home/tema/mdMain/RUT/35/calc/rdf_surf-OW","w")
f_rdf_H = open("/home/tema/mdMain/RUT/35/calc/rdf_surf-HW","w")
for ibin in range(nbins):
    xmax = ax + (ibin+1)*dx
    f_rdf_O.write("{0:10.6f}{1:10.6f}\n".format(xmax, rdf_O[ibin]))
    f_rdf_H.write("{0:10.6f}{1:10.6f}\n".format(xmax, rdf_H[ibin]))
f_rdf_O.close()
f_rdf_H.close()
