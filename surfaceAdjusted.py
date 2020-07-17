'''
Created on Aug 22, 2016

@author: yiqi7710
'''

import numpy as np
import time
import csv
import collections
import math
from scipy.spatial.distance import cdist
import scipy.interpolate as sc
import scipy.spatial as sp
import itertools
import os
import misc
import neighborhood
from metpy.interpolate import interpolate_to_grid

def distance(transect, method, res, path):
    """
    Compute the distance along one transect
    
    Args: 
        - transect (int): index for the transect (e.g., 0 for transect 1)
        - method (string): which method to use to calculate distance
        - dem (tuple): a tuple containing elev. array, geotransform, & cell size
        - path (string): path to the data for the specific study area
    
    Returns:
        - distance along the transect (double)
    """
    # record starting time
    start_time = time.time()

    headerfile= open(path + 'buffer/dem'+str(res)+'m_header.csv', 'rt')
    header = list(csv.reader(headerfile, delimiter=','))[0]

    cellsize=(float(header[2]),-float(header[2]))
    top_left_cor=(float(header[0]),float(header[1]))
    rows=int(header[3])
    cols=int(header[4])
    
    with open(path + 'buffer/'+str(res)+'m/tble_buff_'+str(transect)+'.csv', 'rt') as rasterfile:
        dem=list(csv.reader(rasterfile, delimiter=','))
        dem=np.array(dem[1:len(dem)]).astype(float)

    start_point, end_point = get_start_end_points(path, transect)
    
    if method == "p2p": # pixel to pixel distance
        pnts, elev = p2p_xyz(start_point, end_point, top_left_cor, cellsize, dem)
    else: # must be one of clos, wavg, bilin, biqua, biqub, tin, or nn
        pnts = make_pnts(start_point, end_point, path, cellsize)
        n_pts = pnts.shape[0]
        elev = np.zeros(n_pts)

        if method in ["clos", "wavg", 'biLin', 'biQua', 'biQub']:
            for i in range(n_pts):
                # compute elevation of each sample point
                (x, y) = pnts[i]
                nb_x, nb_y, nb_z = get_nb_vals(i, pnts, dem, top_left_cor, cellsize, rows, cols)
                if method == "clos":
                    elev[i] = nb_z[2, 2]
                elif method == "wavg":
                    elev[i] = neighborhood.WeightedAverage(x, y, nb_x, nb_y, nb_z)
                elif method == "biLin":
                    elev[i] = neighborhood.BiLinear(x, y, nb_x, nb_y, nb_z)
                elif method == "biQua":
                    elev[i] = neighborhood.BiQuadratic(x, y, nb_x, nb_y, nb_z)    
                elif method == "biQub":
                    elev[i] = neighborhood.BiQubic(x, y, nb_x, nb_y, nb_z)
                else:
                    raise ValueError("Unknown method! Must be clos, wavg, biLin, biQua, or biQub")
        else: # must be TIN OR NN
            tn_fname = 'vertices_tin' + str(int(cellsize[0])) + '.csv'
            tn_path=path + tn_fname
            tn_path=tn_path.replace('/simulation','')
            vertices = np.genfromtxt(tn_path, delimiter = ',')
            #Calculate distance between each pair of TIN points (memory consuming)
            #dists = cdist(pnts, vertices[:,0:2], 'euclidean')
            
            #Create a Delaunay triangulation using sample points and TIN vertices.
            tri = sp.Delaunay(vertices[:,0:2])
            #Create Delaunay Triangulation using the sample point and all vertices
            simplex_indices=tri.find_simplex(pnts)
            for i in range(n_pts):
                (x, y) = pnts[i]
                nbr_ver_indx = tri.simplices[simplex_indices[i]]  
                nbr_ver = vertices[nbr_ver_indx]
                #con_ver = get_con_ver(vertices, dists[i], nbr_ver)
                #make an 2D array of the sample sample points because cdist only accepts two arrays
                pnts2=np.array([pnts[i], pnts[i]])
                #Get the distance between the sample point and all other points.
                dists = cdist(pnts2, vertices[:,0:2], 'euclidean')
                #Get nearest vertices around the sample points
                con_ver = get_con_ver(vertices, dists[0], nbr_ver)   
                if method == "TIN":
                    elev[i] = TIN_z(x, y, con_ver, nbr_ver)
                if method == "NN":
                    elev[i] = NN_z(x, y, con_ver, nbr_ver, cellsize)

    dist = 0
    for i in range(len(pnts) - 1):
        d_incr = np.sqrt((pnts[i + 1][0] - pnts[i][0]) ** 2 + 
                         (pnts[i + 1][1] - pnts[i][1]) ** 2 + 
                         (elev[i + 1] - elev[i]) ** 2)
        dist = dist + d_incr
        
    # calculate time elapse
    elapsed_time = time.time() - start_time
    print('time elapsed: '+ str(elapsed_time))
    return dist, elapsed_time


def NN_z(x, y, con_ver, nbr_ver, cellsize):
    """
    Compute elevation using NN method
    """
    gx, gy, elevNNGrid = interpolate_to_grid(con_ver[:, 0], con_ver[:,1], con_ver[:,2], 
                                     interp_type = "natural_neighbor", 
                                     hres = cellsize[0])
    elev_NN = elevNNGrid[0, 0]
    if not(np.isnan(elev_NN)):
        elev_i = elev_NN
    else:
        print("elev_NN is nan: evaluating else loop")
        d_nbr = np.zeros(3)
        for n in range(0, 3):
            d_nbr[n] = ((x - nbr_ver[n][0])**2 + (y - nbr_ver[n][1])**2)**0.5
        nearest_ver = nbr_ver[d_nbr.argmax(0)]
        elev_i = nearest_ver[2]
    return elev_i


def TIN_z(x, y, con_ver, nbr_ver):
    """
    Compute elevation with TIN method
    """
    elev_TIN = sc.griddata((con_ver[:,0], con_ver[:,1]), con_ver[:,2], (x, y), method='linear')
    if not(np.isnan(elev_TIN)):
        elev_i = elev_TIN
    else:
        print("elev_TIN is nan: evaluating else loop")
        d_nbr = np.zeros(3)
        for n in range(0, 3):
            d_nbr[n] = ((x-nbr_ver[n][0])**2 + (y-nbr_ver[n][1])**2)**0.5
        nearest_ver = nbr_ver[d_nbr.argmax(0)]
        elev_i = nearest_ver[2]
    return elev_i


def get_con_ver(vertices, distances, nbr_ver, nbr_pt_num = 100):
    """
    Get the vertices concatenated with nearby points
    """
    dists_ranked = np.sort(distances)
    #if total vertices are less than nbr_pt_num, use all vertices
    if(vertices.shape[0]<nbr_pt_num):
        nearby_pts=vertices
    #Otherwise, use 100 nearest points for NN and 500 for TIN interpolation
    else:
        nearby_pts_idx = distances < dists_ranked[nbr_pt_num]
        nearby_pts=vertices[nearby_pts_idx]
    
    con_ver = np.concatenate((nbr_ver, nearby_pts))
    # remove duplicated points
    con_ver = misc.unique_rows(con_ver)
    return con_ver



def p2p_xyz(start_point, end_point, top_left_cor, cellsize, dem):
    """
    Computes pixel to pixel x, y, and z coordinates for one transect
    """
    start_cell = (int((start_point[0] - top_left_cor[0]) / cellsize[0]),
                  int((start_point[1] - top_left_cor[1]) / cellsize[1]))
    end_cell = (int((end_point[0] - top_left_cor[0]) / cellsize[0]),
                int((end_point[1] - top_left_cor[1]) / cellsize[1]))
    cells = misc.get_line(start_cell, end_cell)    
    pnts = []
    elev = []
    
    dem_elv = dem[:,1]
    dem_indx = dem[:,2:4]

    for cell in cells:
        x = top_left_cor[0] + cell[0] * cellsize[0] + cellsize[0] / 2
        y = top_left_cor[1] + cell[1] * cellsize[1] + cellsize[1] / 2
        #xy_indx=[str(cell[0]),str(cell[1])]
        z_indx=np.logical_and(np.equal(dem_indx[:,0],cell[0]),np.equal(dem_indx[:,1],cell[1]))
        try:
            z=dem_elv[z_indx][0]
        except (np.sum(z_indx)>1):
            print("Oops!  That was more than one indices in dem matching the query index (in getCellValue)")
        #z_indx = [i for i,j in enumerate(dem_indx) if j == xy_indx]
        z = float(dem_elv[z_indx])
        pnts.append((x, y))
        elev.append(z)
    return pnts, elev

def get_start_end_points(path, transect):
    """
    Generates start and end points for a transect
    """        
    transect_array = np.genfromtxt(path + 'tran_sim_pts.csv', delimiter=",")
    start_point = transect_array[2 * transect, :]
    end_point = transect_array[2 * transect + 1, :]
    
    # force start points to be west of end points
    if start_point[0] > end_point[0]:
        previous_start_point = start_point
        start_point = end_point
        end_point = previous_start_point
    return start_point, end_point


def make_pnts(start_point, end_point, path, cellsize):
    """
    Generate sample points for non pixel-to-pixel distance calculations
    """
    pnts = misc.generateSamples(tran = [start_point, end_point], path = path)
    
    #Only do sampling at DEM with >3m resolution 
    if cellsize[0] >= 5:
        pnts = misc.selectSamplePts([start_point, end_point], pnts, cellsize[0])
    return pnts

def get_nb_vals(i, pnts, dem, top_left_cor, cellsize, rows, cols):
    """
    Computes the x coords, y coords, and  elevation in the 
      vicinity of a point in a 5x5 moving window
    """
    nb_x = np.zeros((5,5)) # this 5 by 5 max would contain the x coordinate of 16 neighbor pixels of a sample point
    nb_y = np.zeros((5,5)) # this 5 by 5 matrix would contain the y coordinate of 16 neighbor pixels of a sample point
    nb_z = np.zeros((5,5))
    # get index and value of cell in DEM containing current point
    (cell_X, cell_Y, cell_Z) = misc.getCellValue(pnts[i], 
                                                 dem, 
                                                 top_left_cor, 
                                                 cellsize)
    #Deal with sample points near boundary of the DEM
    point_within_dem = (cell_X-2) >=0 and (cell_Y-2>=0) and (cell_X+3)<=cols and (cell_Y+3)<=rows
    if point_within_dem:
        nb_z[0:5,0:5] = misc.RasterSubset(dem,(cell_Y-2),(cell_Y+3),(cell_X-2),(cell_X+3))
    else:
        #Get the part of moving window within the DEM domain
        in_data= misc.RasterSubset(dem,max((cell_Y-2),0),min((cell_Y+3),rows),max((cell_X-2),0),min((cell_X+3),cols))
        #in_data=dem["array"][max((cell_Y-2),0):min((cell_Y+3),rows),max((cell_X-2),0):min((cell_X+3),cols)]
        nb_z[max((2-cell_Y),0):min((5-(cell_Y+3-rows)),5),max((2-cell_X),0):min((5-(cell_X+3-cols)),5)]=in_data[0:in_data.shape[0],0:in_data.shape[1]]
        in_data_avg=np.mean(in_data[in_data>-3.4e+10])
        nb_z[nb_z==0]=in_data_avg
        nb_z[nb_z<-3.4e+10]=in_data_avg


    
    # If there is missing data in the neighborhood of the sample point 
    # use neighborhood average to replace the missing value 
    has_missing_data = (nb_z>8848).sum()>0 or (nb_z<-413).sum()>0
    if has_missing_data:
        avgValue=np.mean(nb_z[np.where(np.logical_and(nb_z<8848, nb_z>-413))])
        nb_z[nb_z>8848]=avgValue
        nb_z[nb_z<-413]=avgValue
    
    # Obtain the coordinate of cell centroid of a 5*5 neighborhood around the sample point
    for ii in [0,1,2,3,4]:
        cor_y=ii-2
        dy = (cell_Y+cor_y+0.5) * cellsize[1]
        nb_y[ii,:] = top_left_cor[1] + dy
    for jj in [0,1,2,3,4]:
        cor_x=jj-2
        dx = (cell_X+cor_x+0.5) * cellsize[0]
        nb_x [:,jj] = top_left_cor[0] + dx
    return nb_x, nb_y, nb_z


