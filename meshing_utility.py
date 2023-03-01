import pyvista as pv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import vtk
import scipy.stats
import pygalmesh
from vtk.util import numpy_support as npvtk
import json
import shutil

#save data to PSLG file
def mesh2poly(outline, holes=None, hole_pts=None, output_dir='./', scale=1):
    ## inputs:
    ## outline: outerbox boudnary, N*2 dimension, N means total number of points. 
    ## holes: [hole1, hole2,...] hole dimension M*2, M means total number of points of the hole's outline
    ## hole_pts: [hole_pt1, hole_pt2] dimension 2, xy coordinate of pt inside the hole
    ## output_dir: dir of the pslg file, current dir in default
    ## scale: scaling factor, default 1

    path = os.path.join(output_dir, "PSLG")
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


    if holes == None:
        print("internal flow")
        all_pts = outline
        pts = all_pts*scale
        npts = len(pts)
        # starting to write file
        poly_file = open(os.path.join(path, "mesh.poly"), "w")
        poly_file.write("%d %d %d %d\n" % (npts, 2, 1, 0))
        for i in range(npts):
            poly_file.write("%d %f %f %f\n" % (i+1, pts[i,0], pts[i,1], 1.0))
        
        poly_file.write("%d %d\n" % (npts, 1) )
        id_loop= np.concatenate((np.arange(len(outline))+1,[1]))
        for i in range(len(outline)):
            poly_file.write("%d %d %d %d\n" % (i+1, id_loop[i],id_loop[i+1], 0))
        poly_file.write("%d\n" % 0)
        
    else:
        print("number of holes is:",len(holes))
        all_pts = holes.copy()
        all_pts.insert(0,outline)
        # print(holes)
        pts = np.vstack(all_pts)*scale
        hole_pts*=scale
        npts = len(pts)
        # starting to write file
        poly_file = open(os.path.join(path, "mesh.poly"), "w")
        poly_file.write("%d %d %d %d\n" % (npts, 2, 1, 0))
        for i in range(npts):
            poly_file.write("%d %f %f %f\n" % (i+1, pts[i,0], pts[i,1], 1.0))
        
        poly_file.write("%d %d\n" % (npts, 1) )
        id_loop= np.concatenate((np.arange(len(outline))+1,[1]))
        for i in range(len(outline)):
            poly_file.write("%d %d %d %d\n" % (i+1, id_loop[i],id_loop[i+1], 0))
        
        # writing holes
        offset_id = len(outline)
        for ele in holes:
            id_loop= np.concatenate((np.arange(len(ele))+1,[1]))+offset_id
            for i in range(len(ele)):
                poly_file.write("%d %d %d %d\n" % (i+1+offset_id, id_loop[i],id_loop[i+1], 1))
            offset_id += len(ele)
        
        poly_file.write("%d\n" % len(holes))
        for i in range(len(hole_pts)):
            poly_file.write("%d %f %f\n" % (i+1, hole_pts[i][0], hole_pts[i][1]))