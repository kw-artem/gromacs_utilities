import math
import MDAnalysis
import pylab
import numpy as np
from numpy.linalg import norm
import matplotlib as mlab
#import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from matplotlib.ticker import ScalarFormatter

top = "/home/tema/mdMain/RUT/25/md/_rs.gro"
trj = "/home/tema/mdMain/RUT/25/md/_traj.xtc"
chkFound = False
r1 = 2.7; r2 = 6

def plottingHere():
    return None
def toHistogram(array, ax, bx, Np, fdest):
    ylist = np.zeros(Np); k = 0
    ndata = len(array)
    elist = np.sort(array)
    dx = (bx - ax)/Np
    fout = open(fdest, "w")
    fout.write("{0:10d}\n".format(ndata))
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
def seekRound():
    d = coordPart - pos
    r = norm(d, axis=1)
    index = -1
    global chkFound
    val = [None, None, None]
    chkFound = True
    for k in [0, 1, 2]:
        minV = 1000
        for i in range(0,len(r)):
            if r[i] < minV:
                index = i
                minV = r[i]
        if minV > r1:
            chkFound = False
            break
        r[index] = 1000
        val[k] = coordPart[index]
    return val[0], val[1], val[2]
def oxygensInLayer():
    coordOxy = coordSol[::3]
    indices = []
    val = [None, None, None]
    for point in coordPart:
        minV = 1000
        d = coordOxy - point
        r = norm(d, axis=1)
        for i in range(0,len(r)):
            if r[i] < r1:
                chkExist = True
                for j in range(0,len(indices)):
                    if indices[j] == i:
                        chkExist = False
                        break
                if chkExist: indices.append(i)
    return indices
def count(array, cut):
    c = 0
    l = len(array)
    while c < l and array[c] <= cut:
        c += 1
    return c
def addEdge(edge, array):
    a = edge[0]; b = edge[1]
    n = len(array)
    i = 0; j = 0; k = 0
    while i < n and array[i][0] <= a:
        i += 1
    while i < n and array[i][1] > b:
        i -= 1
    array.insert(i, edge)
    return array
def doEdgeExist(edge, array):
    val = True
    #a = edge[0]; b =edge[1]
    for i in range(len(array)):
        if edge == array[i]:
            val = False
            break
    return val
def doFaceExist(face):
    val = True
    #a = edge[0]; b = edge[1]; c = edge[2]
    for i in range(len(faces)):
        if face == faces[i]:
            val = False
            break
    return val
def addTriVertices(face):
    global vertices
    for point in face:
        doExist = False
        for i in range(len(vertices)):
            if point == vertices[i]:
                doExist = True
                break
        if not doExist: vertices.append(point)
    return None
def addTriEdges(face):
    global edges
    sides = [face[:2], face[::2], face[1:]]
    for i in (0, 1, 2):
        if doEdgeExist(sides[i], edges):
            addEdge(sides[i], edges)
    return None
u = MDAnalysis.Universe("md0.gro","trj.trr")

particle = u.select_atoms("name Ti or name O")
solution = u.select_atoms("resname SOL")
R = 1.5
#N_part = len(particle)
#N_sol  = len(solution)
#N_oxy  = N_sol/3
#print "N_sol = %d, N_part = %d" % (N_sol, N_part)
#print(datetime.now())
#N = particle.n_atoms;
#N1 = N/3; N2 = N - N1
#m1 = .47867; m2 = .159994; ms = m1*N1 + m2*N2
#M1 = m1*np.ones(N1, dtype=np.float32)
#M2 = m2*np.ones(N2, dtype=np.float32)
#cm = np.zeros(3, dtype=np.float32)
#def x
#M  = np.concatenate((M1, M2), axis=0)
vertices = []; tEdges = []; edges = []; faces = []
indices  = []
#coordPart = particle.positions
#for i in (0, 1, 2):
#    cm[i] = sum(M*coordPart[:,i])

coordPart = np.array([[0,0,0], [1,0,0], [1,1,0], [0,1,0],
                      [0,0,1], [1,0,1], [1,1,1], [0,1,1]],
                     dtype = np.float32)
cm = np.array([.5,0,0], dtype = np.float32)

d = norm(coordPart - cm, axis=1)
D = (np.sort(d))[::-1]
N = len(coordPart)
x  = np.zeros((N, 3), dtype=np.float32)

for i in range(N):
    curr = D[i]
    k = 0
    while curr != d[k]:
        k += 1
    d[k] = nan
    indices.append(k)
    x[i] = coordPart[k]

for natoms in range(N,0,-1):
    #Edges treatment 
    for i in range(1, natoms):
        act = x[:i]
        dl = norm(act - x[i], axis=1)
        for j in range(i):
            edge = [j, i]
            if dl[j] < R and doEdgeExist(edge, tEdges):
                addEdge(edge, tEdges)
    #Faces treatment
    e = len(tEdges)
    for i in range(e):
        xver = tEdges[i][1]
        for j in range(i+1, e):
            if xver == tEdges[j][0]:
                yver = tEdges[j][1]
                for k in range(i+1, e):
                    if tEdges[k][0] == i and yver == tEdges[k][1]:
                        face = [i, xver, yver]
                        if doFaceExist(face):
                            faces.append(face)
                            addTriVertices(face)
                            addTriEdges(face)
                            break
    v = len(vertices); e = len(edges); f = len(faces)
    chi = v - e + f
    if chi == 2:
        print "The boundary surface is complete!"
        break
print 'done!'
    