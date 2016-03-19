import math
import MDAnalysis
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay

def addVertex(point):
    global boundaryVertices
    boundaryVertices.append(point)
    return None
def doVertexInSet(point):
    val = False
    for i in range(len(boundaryVertices)):
        if point == boundaryVertices[i]:
            val = True
            break
    return val
def addEdge(edge):
    global boundaryEdges
    boundaryEdges.append(edge)
    return None
def doEdgeInSet(edge):
    val = False
    for i in range(len(boundaryEdges)):
        if edge == boundaryEdges[i]:
            val = True
            break
    return val
def whichIsCell(edge):
    val = False
    a = edge[0]; b = edge[0]
    for cell in boundaryCell:
        if (a in cell) and (b in cell):
            val = cell
            break
    return val
def cellFacetSet(cell):
    array = []
    for i in range(4):
        tmpCell = cell
        tmpCell.pop(i)
        array.append(tmpCell)
    return array
def RestructCell(cell, edge):
    global boundaryFacets
    a = edge[0]; b = edge[0]
    i = 0
    for face in boundaryFacets:
        if a in face and b in face:
            boundaryFacets.remove(face)
            i += 1
            if i == 2: break
    facesCell = cellFacetSet(cell)
    for face in facesCell:
        if a in face and b in face:
            facesCell.remove(face)
    boundaryFacets.extend(facesCell)
    return None
DIM  = 3
maxR = 2
#boundaryFacets = np.array([[-1,-1,-1]], dtype = np.int32)
#boundaryEdges  = np.array([[-1,-1]]   , dtype = np.int32)
#boundaryVertices = []
#boundaryFacets   = []
#boundaryEdges    = []
#boundaryCell     = []
#points = np.array([[0,0,0],[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1]])
#points = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]], dtype = np.float32)
#points = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[-1,1,1]])
# Triangulate parameter space to determine the triangles
#tri = mtri.Triangulation(u, v)
u = MDAnalysis.Universe("md0.gro","trj.trr")

particle = u.select_atoms("name Ti or name O")
solution = u.select_atoms("resname SOL")

N_part = len(particle)
N_sol  = len(solution)
N_oxy  = N_sol/3

outPoints = []
points   = particle.positions
init     = particle.positions
iteration = 0
while iteration < 3:
    iteration += 1
    boundaryVertices = []
    boundaryFacets   = []
    boundaryEdges    = []
    boundaryCell     = []
    triCells = Delaunay(points, False)
    allCells = triCells.simplices; allCellsl = allCells.tolist()
    ncells   = triCells.nsimplex
    #for vert in tri.triangles:
    print allCells, ncells
    
    neighborCells  = triCells.neighbors
    
    #Searching for facets at the boundary of the triangulation
    l = -1
    for i in range(ncells):
        atBoundary = False
        for k in range(DIM + 1):
            if neighborCells[i,k] == -1:
                l += 1
                atBoundary = True
                tmpRow = np.delete(allCells[i], k).tolist()
                boundaryFacets.append(tmpRow)
                #boundaryFacets = np.row_stack((boundaryFacets, tmpRow))
        if atBoundary: boundaryCell.append(allCells[i])
    nbfacets = len(boundaryFacets)# - 1
    #print boundaryFacets, nbfacets
    print nbfacets
    
    #Forming a set of edges from some lies at the boundary
    for i in range(nbfacets):
        face  = boundaryFacets[i]
        sides = [face[:2],face[1:],face[::2]]
        for j in range(3):
            if not doEdgeInSet(sides[j]): addEdge(sides[j])
    nbedges = len(boundaryEdges)
    print boundaryEdges, nbedges
    
    #Forming a set of points which are surface vertices
    fout = open("_idsSurf.txt", "w")#
    for i in range(nbedges):
        segment  = boundaryEdges[i]
        for j in range(2):
            if not doVertexInSet(segment[j]): 
                addVertex(segment[j])
                fout.write(" {0}".format(segment[j]))#
    fout.close()#
    nbvertices = len(boundaryVertices)
    print boundaryVertices, nbvertices
    
    coordbEdges = np.zeros((nbedges, DIM-1), dtype = np.float32)
    #Checking whether each boundary edge satisfies the distance condition
    tmpArr_chk = np.ones((nbedges), dtype = np.int32)
    for i in range(nbedges):
        segment = points[boundaryEdges[i]]
        length  = norm(segment[0] - segment[1])
        if length > maxR:
            tmpArr_chk[i] = 0
    
    boundaryVertices.sort()
    for i in boundaryVertices:
        point = points[i].tolist()
        for j in range(len(init)):
            if point == init[j].tolist():
                outPoints.append(j)
                break
    outPoints.sort()
#    if iteration < 2:
#        outPoints = boundaryVertices
#        boundaryVertices.sort()
#        points = np.delete(points, boundaryVertices, 0)
    points = np.delete(points, boundaryVertices, 0)



fout = open("_idsSurf.txt", "w")
for i in outPoints:
    fout.write(" {0}".format(i))
fout.close()
'''
#Removing bad edges
for i in range(len(tmpArr_chk)):
    if tmpArr_chk[i] == 0:
        edge = boundaryEdges.pop(i)
        cell = whichIsCell(edge)
        RestructCell(cell, edge)
nbfacets = len(boundaryFacets)# - 1

print "++++++++++++++++++\nBefore"
print boundaryFacets, nbfacets

#Forming a set of edges from some lies at the boundary
for i in range(nbfacets):
    face  = boundaryFacets[i]
    sides = [face[:2],face[1:],face[::2]]
    for j in range(3):
        if not doEdgeInSet(sides[j]): addEdge(sides[j])
nbedges = len(boundaryEdges)
print boundaryEdges, nbedges

#Forming a set of points which are surface vertices
fout = open("_idsSurf.txt", "w")#
for i in range(nbedges):
    segment  = boundaryEdges[i]
    for j in range(2):
        if not doVertexInSet(segment[j]): 
            addVertex(segment[j])
            fout.write(" {0}".format(segment[j]))#
fout.close()#
nbvertices = len(boundaryVertices)
print boundaryVertices, nbvertices

coordbEdges = np.zeros((nbedges, DIM-1), dtype = np.float32)
#Checking whether each boundary edge satisfies the distance condition
tmpArr_chk = np.ones((nbedges), dtype = np.int32)
for i in range(nbedges):
    segment = points[boundaryEdges[i]]
    length  = norm(segment[0] - segment[1])
    if length > maxR:
        tmpArr_chk[i] = 0
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1, projection='3d')

# The triangles in parameter space determine which x, y, z points are
# connected by an edge
#ax.plot_trisurf(x, y, z, triangles=tri.triangles, cmap=plt.cm.Spectral)
#ax.plot_trisurf(x, y, z, triangles=tri.simplices, cmap=plt.cm.Spectral)
#plt.show()
'''