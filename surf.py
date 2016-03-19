import math
import MDAnalysis
import numpy as np
from numpy.linalg import norm
from numpy.random import random_sample
from datetime import datetime, timedelta
import matplotlib.pyplot as pl
import sys

def doIndexExist(element, array):
    val = False
    for z in array:
        if element == z:
            val = True
            break
    return val
def sortMax(array, indices):
    maxval = array[0]
    index  = indices[0]
    for i in range(len(array)):
        if array[i] > maxval:
            maxval = array[i]
            index  = indices[i]
    return index
def sliceSurfaceAtoms(coordPart, r = 1.4, D = 8, n = 3, p = 1):
    #R=10,r=1.4,D=8,n=2,p=1
    N = len(coordPart)
    N1 = N/3; N2 = N - N1
    #
    #m1 = .47867; m2 = .159994; 
    #ms = m1*N1 + m2*N2
    #
    #M1 = m1*np.ones(N1, dtype=np.float32)
    #M2 = m2*np.ones(N2, dtype=np.float32)
    #cm = np.zeros(3, dtype=np.float32)
    #M  = np.concatenate((M1, M2), axis=0)
    #
    coordPart = particle.positions
    cm = coordPart.mean(axis=0)
    #cm = M.dot(coordPart)/ms
    #print cm
    dev = 0
    nloop = 0
    print(datetime.now())
    while nloop < n:
        x0 = cm + dev
        dx = coordPart - x0
        dl = norm(dx, axis=1)
        
        if nloop == 0:
            k = 0
            indices = []
            surface = []
            for rl in dl:
                if rl > .8*D/2:
                    indices.append(k)
                k += 1
            ncurr   = len(indices)
            points  = coordPart[indices]
            dists   = dl[indices]
            weigths = np.zeros(N, dtype = np.int32)
        
        
        for inda in range(N):
            an = dx[inda]
            atomsInRay = []
            for indb in indices:
                bn = dx[indb]
                cos = np.dot(an, bn)/dl[inda]/dl[indb]
                MIN = 1/sqrt((1 + tan(r/dl[indb])**2))
                if cos >= MIN:
                    atomsInRay.append(indb)
            distInRay  = dl[atomsInRay]
            farthest   = sortMax(distInRay, atomsInRay)
            weigths[farthest] += 1 
    
        nloop += 1
        if n > 1:
            dev = .8*D/2*random_sample((3,))
    print(datetime.now())
    
    for i in range(N):
        if weigths[i] > p:
            surface.append(i)
    return surface

if __name__ == "__main__":
    _indfile = True
    try:
        topology   = sys.argv[1]
        trajectory = sys.argv[2]
    except IndexError:
        topology   = 'md0.gro'
        trajectory = 'trj.trr'
        indicesOut = 'idsSurf.txt'
    try:
        indicesOut = sys.argv[3]
    except IndexError:
        if not indicesOut: _indfile = False
                
    particle = u.select_atoms("name Ti or name O")
    surface  = sliceSurfaceAtoms(particle.positions)
    if _indfile:
        fout = open(indicesOut, "w")
        for ind in surface:
            fout.write(" {0}".format(ind))
        fout.close()
