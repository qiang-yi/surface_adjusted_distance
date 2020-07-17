'''
Created on Sep 9, 2016

@author: yiqi7710
'''

import numpy as np
import local_interp

def WeightedAverage(x,y,rasterBlock_x,rasterBlock_y, rasterBlock_elev):
    
    xCoor = np.array([rasterBlock_x [2,1], rasterBlock_x[1,2], rasterBlock_x[2,3], rasterBlock_x[3,2], rasterBlock_x[1,1],
             rasterBlock_x[1,3], rasterBlock_x[2,2], rasterBlock_x[3,1], rasterBlock_x[3,3]])
    yCoor = np.array([rasterBlock_y [2,1], rasterBlock_y[1,2], rasterBlock_y[2,3], rasterBlock_y[3,2], rasterBlock_y[1,1],
             rasterBlock_y[1,3], rasterBlock_y[2,2], rasterBlock_y[3,1], rasterBlock_y[3,3]])
    elev = np.array([rasterBlock_elev [2,1], rasterBlock_elev[1,2], rasterBlock_elev[2,3], rasterBlock_elev[3,2], rasterBlock_elev[1,1],
             rasterBlock_elev[1,3], rasterBlock_elev[2,2], rasterBlock_elev[3,1], rasterBlock_elev[3,3]])
    elevWeiAvr=local_interp.IDW(x,y, xCoor, yCoor, elev, 2)
    return elevWeiAvr

def BiLinear(x,y,rasterBlock_x,rasterBlock_y, rasterBlock_elev):
    
    xCoor = np.array([rasterBlock_x [2,1], rasterBlock_x[1,2], rasterBlock_x[2,3], rasterBlock_x[3,2],rasterBlock_x[2,2]])
    yCoor = np.array([rasterBlock_y [2,1], rasterBlock_y[1,2], rasterBlock_y[2,3], rasterBlock_y[3,2],rasterBlock_y[2,2]])
    elev = np.array([rasterBlock_elev [2,1], rasterBlock_elev[1,2], rasterBlock_elev[2,3], rasterBlock_elev[3,2],rasterBlock_elev[2,2]])
    elevBiLinear=local_interp.polyval2d(x,y,local_interp.polyfit2d(xCoor, yCoor, elev, 1))
    return elevBiLinear

def BiQuadratic(x,y,rasterBlock_x,rasterBlock_y, rasterBlock_elev):
    
    xCoor = np.array([rasterBlock_x [2,1], rasterBlock_x[1,2], rasterBlock_x[2,3], rasterBlock_x[3,2], rasterBlock_x[1,1],
             rasterBlock_x[1,3], rasterBlock_x[2,2], rasterBlock_x[3,1], rasterBlock_x[3,3]])
    yCoor = np.array([rasterBlock_y [2,1], rasterBlock_y[1,2], rasterBlock_y[2,3], rasterBlock_y[3,2], rasterBlock_y[1,1],
             rasterBlock_y[1,3], rasterBlock_y[2,2], rasterBlock_y[3,1], rasterBlock_y[3,3]])
    elev = np.array([rasterBlock_elev [2,1], rasterBlock_elev[1,2], rasterBlock_elev[2,3], rasterBlock_elev[3,2], rasterBlock_elev[1,1],
             rasterBlock_elev[1,3], rasterBlock_elev[2,2], rasterBlock_elev[3,1], rasterBlock_elev[3,3]])
    elevBiQuadratic=local_interp.polyval2d(x,y,local_interp.polyfit2d(xCoor, yCoor, elev, 2))
    return elevBiQuadratic

def BiQubic(x,y,rasterBlock_x,rasterBlock_y, rasterBlock_elev):
    
    xCoor = np.array([rasterBlock_x [2,1], rasterBlock_x[1,2], rasterBlock_x[2,3], rasterBlock_x[3,2], rasterBlock_x[1,1],
            rasterBlock_x[1,3], rasterBlock_x[3,1], rasterBlock_x[3,3], rasterBlock_x[0,0], rasterBlock_x[0,2], rasterBlock_x[0,4],
            rasterBlock_x[2,0], rasterBlock_x[2,2], rasterBlock_x[2,4], rasterBlock_x[4,0], rasterBlock_x[4,2], rasterBlock_x[4,4]])
    yCoor = np.array([rasterBlock_y [2,1], rasterBlock_y[1,2], rasterBlock_y[2,3], rasterBlock_y[3,2], rasterBlock_y[1,1],
            rasterBlock_y[1,3], rasterBlock_y[3,1], rasterBlock_y[3,3], rasterBlock_y[0,0], rasterBlock_y[0,2], rasterBlock_y[0,4],
            rasterBlock_y[2,0], rasterBlock_y[2,2], rasterBlock_y[2,4], rasterBlock_y[4,0], rasterBlock_y[4,2], rasterBlock_y[4,4]])
    elev = np.array([rasterBlock_elev [2,1], rasterBlock_elev[1,2], rasterBlock_elev[2,3], rasterBlock_elev[3,2], rasterBlock_elev[1,1],
            rasterBlock_elev[1,3], rasterBlock_elev[3,1], rasterBlock_elev[3,3], rasterBlock_elev[0,0], rasterBlock_elev[0,2], rasterBlock_elev[0,4],
            rasterBlock_elev[2,0], rasterBlock_elev[2,2], rasterBlock_elev[2,4], rasterBlock_elev[4,0], rasterBlock_elev[4,2], rasterBlock_elev[4,4]])
    elevBiQubic=local_interp.polyval2d(x,y,local_interp.polyfit2d(xCoor, yCoor, elev, 3))
    return elevBiQubic

