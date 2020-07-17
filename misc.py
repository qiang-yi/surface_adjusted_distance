'''
Created on Aug 25, 2016




@author: yiqi7710
'''
# this finction creates ragular sample points along the input polyline based on the sample size
from __future__ import division
import math
import collections
import scipy.spatial
import numpy as np
import time


# def getNearestNeighbors(pt,top_left_cor,cellsize):
#     ctr_cell=(int((pt[0]-top_left_cor[0])/cellsize[0]),int((pt[1]-top_left_cor[1])/cellsize[1]))
#     toLeft=pt[0]-ctr_cell[0]*cellsize[0]
#     toTop=pt[1]-ctr_cell[1]*cellsize[1]
#     if abs(toLeft)<0.5*abs(cellsize[0]) and abs(toTop)<0.5*abs(cellsize[1]):
#         celltl=((ctr_cell[0]-1),ctr_cell[1]+1)
#         celltr=((ctr_cell[0]),ctr_cell[1]+1)
#         cellbl=((ctr_cell[0]-1),ctr_cell[1])
#         celltr=ctr_cell
#     elif abs(toLeft)>0.5*abs(cellsize[0]) and abs(toTop)<0.5*abs(cellsize[1]):
#         celltl=((ctr_cell[0]),ctr_cell[1]+1)
#         celltr=((ctr_cell[0]+1),ctr_cell[1]+1)
#         cellbl=ctr_cell
#         celltr=((ctr_cell[0]+1),ctr_cell[1])

def getCellValue(point,dem,top_left_cor,cellsize):
    X=int((point[0]-top_left_cor[0])/cellsize[0])
    Y=int((point[1]-top_left_cor[1])/cellsize[1])
    
    dem_elv = dem[:,1]
    dem_indx = dem[:,2:4]
    z_indx=np.logical_and(np.equal(X,dem_indx[:,0]),np.equal(Y,dem_indx[:,1]))
#    dem_indx.index(np.array([X,Y]))
    try:
        Z=dem_elv[z_indx][0]
    except (np.sum(z_indx)>1):
        print("Oops!  That was more than one indices in dem matching the query index (in getCellValue)")
    
    return (X,Y,Z)

def RasterSubset(dem,bottom,top,left,right):
    X_inbox=np.logical_and(np.greater_equal(dem[:,2],left),np.less(dem[:,2],right))
    Y_inbox=np.logical_and(np.greater_equal(dem[:,3],bottom),np.less(dem[:,3],top))
    XY_inbox=np.logical_and(X_inbox,Y_inbox)
    dem_inbox=dem[XY_inbox]
    
    indx_list=np.array([[x,y] for y in range(bottom,top) for x in range(left,right)])

    if (len(dem_inbox)==len(indx_list)):
        shape_dif_sum=np.sum(abs(np.array(indx_list.shape[0])-np.array(dem_inbox.shape[0])))
        indx_dif_sum=np.sum(abs(np.array(indx_list)-np.array(dem_inbox[:,2:4])))
        
        if (shape_dif_sum>0):
            print("Error, shape of dem_inbox and moving window does not match (in RasterSubset)!")
            exit()
        #if (indx_dif_sum>0):
            #print("Error, XY indices of dem_inbox and moving window does not match (in RasterSubset)!")
    
        z_list=dem_inbox[:,1]
        z_chuck=np.reshape(z_list,((top-bottom),(right-left)))
    else:
        #print('There is missing data in the moving window. bottom:'+str(bottom)+ ' top: '+str(top)+' left: '+str(left)+ ' right: '+str(right))
        dem_indx=[str(i)+'_'+str(j) for i,j in list(dem[:,2:4].astype(int))]
        z_chuck=np.zeros(((top-bottom),(right-left)))
        for y in range(bottom,top):
            for x in range(left,right):
                xy_indx=str(x)+'_'+str(y)
                if xy_indx in dem_indx:
                    indx=dem_indx.index(xy_indx)
                    z_chuck[(y-bottom),(x-left)]=dem[indx,1]
                else:
                    z_chuck[(y-bottom),(x-left)]=np.nan
        avg=np.nanmean(z_chuck)
        z_chuck[np.isnan(z_chuck)]=avg
    #print "processing time for getting a raster truck is: "+str(time.time()-t)

    return z_chuck

def generateSamples(tran,path):
    sampleSize=3
    #tran_FID=tran['properties']['Id']
    start_pt=tran[0]
    end_pt=tran[1]
    
    x_len=end_pt[0]-start_pt[0]
    y_len=end_pt[1]-start_pt[1]
    tran_len=math.sqrt(math.pow((end_pt[1]-start_pt[1]),2)+math.pow((end_pt[0]-start_pt[0]),2))
    cos=x_len/tran_len
    sin=y_len/tran_len
    pt=start_pt
    pts=[]
    cur_len=0
    #Add the first point
    pts.append(start_pt)
    while True:
        x_inc=sampleSize*cos
        y_inc=sampleSize*sin
        pt=((pt[0]+x_inc),(pt[1]+y_inc))
        cur_len=cur_len+sampleSize
        #When the current length of sample pt string is longer than the transect, jump out the loop
        if cur_len>tran_len:
            break
        else:
            pts.append(pt)
    #Finally, add the last sample point
    pts.append(end_pt)
    
    #Load 3m DEM and stdv of 3m DEM
#     ds1=gdal.Open(path+"dem3m")
#     dem3m=ds1.GetRasterBand(1).ReadAsArray()
#     ds2=gdal.Open(path+"dem3m_st")
#     dem3m_st=ds2.GetRasterBand(1).ReadAsArray()
    
    #Get the resolutions of DEM list to name the table columns
#     gt = ds1.GetGeoTransform()
#     cellsize=(gt[1],gt[5])
#     top_left_cor=(gt[0],gt[3])
#     dem_ls=[]
#     dem_st_ls=[]
#     
#     for pt in pts:
#         dem=misc.getCellValue(pt,dem3m,top_left_cor,cellsize)[2]
#         dem_ls.append(dem)
#         dem_st=misc.getCellValue(pt,dem3m_st,top_left_cor,cellsize)[2]
#         dem_st_ls.append(dem_st)
#     
    #Zip point list into [(x_cor,y_cor),dem,dem_st]
#     pts=zip(pts,dem_ls,dem_st_ls)
#     
#     pts_fn='s_pts_tran'+str(tran_FID)+'_size'+str(int(sampleSize))+'.shp'
    #pts_fn='aaaa.shp'
#     print('output file: '+pts_fn)
    
    #Create a point shapefile to store the sample point
    #Copy the schema from transects shapefile
#     new_crs={'init': u'epsg:26913'}
#     new_driver='ESRI Shapefile'
#     new_schema={'geometry': 'Point', 'properties': collections.OrderedDict([(u'Id', 'int'), (u'x_cor', 'float'), (u'y_cor', 'float'), (u'dem', 'float'), (u'dem_st', 'float')])}
#     points=fiona.open(path+pts_fn,'w',driver=new_driver,crs=new_crs,schema=new_schema)        
#     for i in range(0,len(pts)):
#         rec={'geometry': {'type': 'Point', 'coordinates': (pts[i][0][0],pts[i][0][1])}, 'type': 'Feature', 'id': '0', 'properties': collections.OrderedDict([(u'Id', i), (u'x_cor', pts[i][0][0]), (u'y_cor', pts[i][0][1]), (u'dem', pts[i][1]), (u'dem_st', pts[i][2])])}
#         points.write(rec)
#     
#     points.close()
    pts=np.asarray(pts)
    return pts


def selectSamplePts(tran,samplePts,resolution):
    start_pt=tran[0]
    end_pt=tran[1]
    tran_len=((end_pt[0]-start_pt[0])**2+(end_pt[1]-start_pt[1])**2)**0.5
    num_pts=int(tran_len/resolution)*3
    pool_size=samplePts.shape[0]
    num_pts=min(pool_size,num_pts)
    rand_indx=np.random.choice(pool_size, num_pts,replace=False)
    rand_indx=np.sort(rand_indx)
    pts=samplePts[rand_indx]
    return pts

def randomSampling(tran,resolution,path):
    #tran_FID=tran['properties']['Id']
    start_pt=tran[0]
    end_pt=tran[1]
    
    x_len=end_pt[0]-start_pt[0]
    y_len=end_pt[1]-start_pt[1]
    tran_len=((end_pt[0]-start_pt[0])**2+(end_pt[1]-start_pt[1])**2)**0.5
    slope=y_len/x_len   

    pt_per_cell=3
    num_pts=tran_len/resolution*pt_per_cell
    x_incr=x_len*np.random.random(num_pts)
    x_incr=np.sort(x_incr)
    Xs=start_pt[0]+x_incr
    Ys=start_pt[1]+slope*x_incr
    pts=np.column_stack((Xs,Ys))
        
    #Load 3m DEM and stdv of 3m DEM
#     ds1=gdal.Open(path+"dem3m")
#     dem3m=ds1.GetRasterBand(1).ReadAsArray()
#     ds2=gdal.Open(path+"dem3m_st")
#     dem3m_st=ds2.GetRasterBand(1).ReadAsArray()
    
    #Get the resolutions of DEM list to name the table columns
#     gt = ds1.GetGeoTransform()
#     cellsize=(gt[1],gt[5])
#     top_left_cor=(gt[0],gt[3])
#     dem_ls=[]
#     dem_st_ls=[]
#     
#     for pt in pts:
#         dem=misc.getCellValue(pt,dem3m,top_left_cor,cellsize)[2]
#         dem_ls.append(dem)
#         dem_st=misc.getCellValue(pt,dem3m_st,top_left_cor,cellsize)[2]
#         dem_st_ls.append(dem_st)
#     
    #Zip point list into [(x_cor,y_cor),dem,dem_st]
#     pts=zip(pts,dem_ls,dem_st_ls)
#     
    #pts_fn='s_pts_tran'+str(tran_FID)+'_size'+str(int(resolution))+'.shp'
    #print('output file: '+pts_fn)
    
    #Create a point shapefile to store the sample point
    #Copy the schema from transects shapefile
#     new_crs={'init': u'epsg:26913'}
#     new_driver='ESRI Shapefile'
#     new_schema={'geometry': 'Point', 'properties': collections.OrderedDict([(u'Id', 'int'), (u'x_cor', 'float'), (u'y_cor', 'float')])}
#     points=fiona.open(r"C:/Users/yiqi7710/work_CU/data/Modified_2/Nebraska/sample_points/"+pts_fn,'w',driver=new_driver,crs=new_crs,schema=new_schema)
#     
#     for i in range(0,num_pts):
#         rec={'geometry': {'type': 'Point', 'coordinates': (pts[i][0],pts[i][1])}, 'type': 'Feature', 'id': '0', 'properties': collections.OrderedDict([(u'Id', i), (u'x_cor', pts[i][0]), (u'y_cor', pts[i][1])])}
#         points.write(rec)
#      
#     points.close()
    return pts

def get_line(start, end):
    """Bresenham's Line Algorithm
    Produces a list of tuples from start and end
 
    >>> points1 = get_line((0, 0), (3, 4))
    >>> points2 = get_line((3, 4), (0, 0))
    >>> assert(set(points1) == set(points2))
    >>> print points1
    [(0, 0), (1, 1), (1, 2), (2, 3), (3, 4)]
    >>> print points2
    [(3, 4), (2, 3), (1, 2), (1, 1), (0, 0)]
    """
    # Setup initial conditions
    x1, y1 = start
    x2, y2 = end
    dx = x2 - x1
    dy = y2 - y1
 
    # Determine how steep the line is
    is_steep = abs(dy) > abs(dx)
 
    # Rotate line
    if is_steep:
        x1, y1 = y1, x1
        x2, y2 = y2, x2
 
    # Swap start and end points if necessary and store swap state
    swapped = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        swapped = True
 
    # Recalculate differentials
    dx = x2 - x1
    dy = y2 - y1
 
    # Calculate error
    error = int(dx / 2.0)
    ystep = 1 if y1 < y2 else -1
 
    # Iterate over bounding box generating points between start and end
    y = y1
    points = []
    for x in range(x1, x2 + 1):
        coord = (y, x) if is_steep else (x, y)
        points.append(coord)
        error -= abs(dy)
        if error < 0:
            y += ystep
            error += dx
 
    # Reverse the list if the coordinates were swapped
    if swapped:
        points.reverse()
    return points


def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))