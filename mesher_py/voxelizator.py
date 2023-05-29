#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 11:41:27 2022

@author: anonym
"""
# import convert
from stl.mesh import Mesh
import numpy as np
from .voxelize_internal import voxel_intern
from numba import jit

@jit(nopython=True, cache=True, fastmath=True)
def solve_overlapping(n_cells_x,n_cells_y,n_cells_z,num_materials,id_mat_keep,output_meshing):

    for c1 in range(n_cells_x):
        for c2 in range(n_cells_y):
            for c3 in range(n_cells_z):
                for k in range(num_materials):
                    if output_meshing[k,c1, c2, c3] == True:
                        for k2 in range(num_materials):
                            if output_meshing[k2,c1, c2, c3] == True and k!=k2:
                                if k in id_mat_keep:
                                    output_meshing[k2, c1, c2, c3] = False
                                else:
                                    output_meshing[k, c1, c2, c3] = False

@jit(nopython=True, cache=True, fastmath=True)
def merge_the_3_grids(voxcountX, voxcountY, voxcountZ,gridOUTPUT1,gridOUTPUT2,gridOUTPUT3,gridOUTPUT):
    for c1 in range(voxcountX):
        for c2 in range(voxcountY):
            for c3 in range(voxcountZ):
                cont_true = 0
                if gridOUTPUT1[c2, c3, c1] == True:
                    cont_true = cont_true + 1
                if gridOUTPUT2[c3, c1, c2] == True:
                    cont_true = cont_true + 1
                if gridOUTPUT3[c1, c2, c3] == True:
                    cont_true = cont_true + 1
                if cont_true > 1:
                    gridOUTPUT[c1, c2, c3] = True

def voxelize(cells_on_x:int,cells_on_y:int,cells_on_z:int,meshXYZ:Mesh, geometry_desc:dict):
    
    """
    geometry_dec è il vecchio minmax, lo passiamo come parametro, evitare variabili globali
    """
    assert type(cells_on_x)==int
    assert type(cells_on_y)==int
    assert type(cells_on_z)==int
    assert type(geometry_desc)==dict

    assert len(geometry_desc)==6
    
    # raydirection    = 'xyz'
    
    assert type(meshXYZ)==Mesh
    
    meshXmin = geometry_desc['meshXmin']
    meshXmax = geometry_desc['meshXmax']
    meshYmin = geometry_desc['meshYmin']
    meshYmax = geometry_desc['meshYmax']
    meshZmin = geometry_desc['meshZmin']
    meshZmax = geometry_desc['meshZmax']


    if cells_on_x==1:
        #If gridX is a single integer (rather than a vector) and is equal to 1
        gridCOx   = (meshXmin+meshXmax)/2
    else: 
        #If gridX is a single integer (rather than a vector) then automatically create the list of x coordinates
        voxwidth  = (meshXmax-meshXmin)/(cells_on_x+1/2)
        gridCOx   = np.arange(meshXmin+voxwidth/2, meshXmax-voxwidth/2, voxwidth)
    
    if cells_on_y==1:
        #If gridX is a single integer (rather than a vector) and is equal to 1
        gridCOy   = (meshYmin+meshYmax)/2
    else: 
        #If gridX is a single integer (rather than a vector) then automatically create the list of x coordinates
        voxwidth  = (meshYmax-meshYmin)/(cells_on_y+1/2)
        gridCOy   = np.arange(meshYmin+voxwidth/2, meshYmax-voxwidth/2, voxwidth)
    
    if cells_on_z==1:
        #If gridX is a single integer (rather than a vector) and is equal to 1
        gridCOz   = (meshZmin+meshZmax)/2
    else: 
        #If gridX is a single integer (rather than a vector) then automatically create the list of x coordinates
        voxwidth  = (meshZmax-meshZmin)/(cells_on_z+1/2)
        gridCOz   = np.arange(meshZmin+voxwidth/2, meshZmax-voxwidth/2, voxwidth)

    assert min(gridCOx)>meshXmin
    # Check that the output grid is large enough to cover the mesh:
    var_check=0
    if (min(gridCOx)>meshXmin or max(gridCOx)<meshXmax):
        var_check = 1
        gridcheckX = 0
        if min(gridCOx)>meshXmin:
            gridCOx=np.hstack((meshXmin, gridCOx))
            gridcheckX = gridcheckX+1
        if max(gridCOx)<meshXmax:
            gridCOx = np.hstack((gridCOx,meshXmax))
            gridcheckX = gridcheckX+2
    else:
        if (min(gridCOy)>meshYmin or max(gridCOy)<meshYmax):
            var_check = 2
            gridcheckY = 0
            if min(gridCOy)>meshYmin:
                gridCOy = np.hstack((meshYmin, gridCOy))
                gridcheckY = gridcheckY+1
            if max(gridCOy)<meshYmax:
                gridCOy = np.hstack((gridCOy, meshYmax))
                gridcheckY = gridcheckY+2
        else:
            if (min(gridCOz)>meshZmin or max(gridCOz)<meshZmax):
                var_check = 3
                gridcheckZ = 0
                if min(gridCOz)>meshZmin:
                    gridCOz = np.hstack((meshZmin, gridCOz))
                    gridcheckZ = gridcheckZ+1
                if max(gridCOz)<meshZmax:
                    gridCOz = np.hstack((gridCOz, meshZmax))
                    gridcheckZ = gridcheckZ+2

    #print(gridCOx)
    #print(gridCOy)
    #print(gridCOz)

    # %======================================================
    # % VOXELISE USING THE USER DEFINED RAY DIRECTION(S)
    # %======================================================
    
    # Count the number of voxels in each direction:
    voxcountX = len(gridCOx)
    voxcountY = len(gridCOy)
    voxcountZ = len(gridCOz)

    gridOUTPUT1 = voxel_intern(gridCOy, gridCOz, gridCOx, meshXYZ.v0, meshXYZ.v1, meshXYZ.v2, geometry_desc, 0)
    gridOUTPUT2 = voxel_intern(gridCOz, gridCOx, gridCOy, meshXYZ.v0, meshXYZ.v1, meshXYZ.v2, geometry_desc, 1)
    gridOUTPUT3 = voxel_intern(gridCOx, gridCOy, gridCOz, meshXYZ.v0, meshXYZ.v1, meshXYZ.v2, geometry_desc, 2)

    #from scipy.io import savemat
    #mdic = {"gridOUTPUT1": gridOUTPUT1, "gridOUTPUT2": gridOUTPUT2, "gridOUTPUT3": gridOUTPUT3}
    #savemat("python_grids.mat", mdic)

    gridOUTPUT = np.full((voxcountX, voxcountY, voxcountZ), False, dtype=bool)
    merge_the_3_grids(voxcountX, voxcountY, voxcountZ, gridOUTPUT1, gridOUTPUT2, gridOUTPUT3,gridOUTPUT)

    # %======================================================
    # % RETURN THE OUTPUT GRID TO THE SIZE REQUIRED BY THE USER (IF IT WAS CHANGED EARLIER)
    # %======================================================
    assert var_check in [1,2,3]
    # match var_check:
    if var_check == 1:
        assert gridcheckX in [1,2,3]

        if gridcheckX == 1:
            gridOUTPUT=gridOUTPUT[1:voxcountX,:,:]
            gridCOx    = gridCOx[1:voxcountX]
        elif gridcheckX == 2:
            gridOUTPUT = gridOUTPUT[0:voxcountX-1, :, :]
            gridCOx = gridCOx[0:voxcountX-1]
        else:
            assert gridcheckX == 3
            gridOUTPUT = gridOUTPUT[1:voxcountX - 1, :, :]
            gridCOx = gridCOx[1:voxcountX - 1]
    elif var_check == 2:
            assert gridcheckY in [1,2,3]

            if gridcheckY == 1:
                gridOUTPUT = gridOUTPUT[:,1:voxcountY, :]
                gridCOy = gridCOy[1:voxcountY]
            elif gridcheckY == 2:
                gridOUTPUT = gridOUTPUT[:, 0:voxcountY-1, :]
                gridCOy = gridCOy[0:voxcountY-1]
            else:
                assert gridcheckY == 3
                gridOUTPUT = gridOUTPUT[:, 1:voxcountY-1, :]
                gridCOy = gridCOy[1:voxcountY-1]
    else:
        assert var_check == 3
        assert gridcheckZ in [1,2,3]
        if gridcheckZ == 1:
            gridOUTPUT = gridOUTPUT[:, :, 1:voxcountZ]
            gridCOz = gridCOz[1:voxcountZ]
        elif gridcheckZ == 2:
            gridOUTPUT = gridOUTPUT[:, :, 0:voxcountZ-1]
            gridCOz = gridCOz[0:voxcountZ-1]
        else:
            assert gridcheckZ == 3
            gridOUTPUT = gridOUTPUT[:, :, 1:voxcountZ-1]
            gridCOz = gridCOz[1:voxcountZ-1]

    #if needed also gridCOx, gridCOy and gridCOz can be given in output
    return gridOUTPUT