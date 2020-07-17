'''
Created on Sep 13, 2016

@author: yiqi7710
'''
import itertools
import numpy as np
import matplotlib.mlab as ml
import math

# This module can be used to fit a polynomial surface to a set of (x,y,z) points
#Reference: http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent


def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z

def IDW(x, y, xCoords, yCoords, elevs, order=2):

    d =  np.zeros(len(xCoords))
    w =  np.zeros(len(xCoords))
    for i in range(len(xCoords)):
        d[i] = math.sqrt((x-xCoords[i])**2 +(y-yCoords[i])**2)
        w[i] = 1.0 / d[i]**order

    z = (np.sum(elevs*w)) / float(np.sum(w))
    return z

# def NaturalNeighbor(vertices,dist_to_ver_i,x,y,cellsize):
#     dists_ranked=np.sort(dist_to_ver_i)
#     nearby_pts = dist_to_ver_i < dists_ranked[100]
#     if(nearby_pts[nearby_pts==True].sum()<100):
#         nearby_pts=vertices
#     elevNNGridX=np.array([x,y+cellsize[0]])
#     elevNNGridY=np.array([x,y+cellsize[0]])
#     
#     elevNNGrid=ml.griddata(vertices[nearby_pts,0], vertices[nearby_pts,1], vertices[nearby_pts,2], elevNNGridX, elevNNGridY, interp ='nn')
#     return elevNNGrid[0][0]

