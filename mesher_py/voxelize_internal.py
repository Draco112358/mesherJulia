# -*- coding: utf-8 -*-
"""
Created on Fri May  6 08:22:21 2022

@author: anonym
"""
import numpy as np
from numpy import linalg as LA
from numba import jit

@jit(nopython=True, cache=True, fastmath=True)
def find_pos_min(vect,N):
    min_value=1e30
    pos=-1
    for cont in range(N):
        if vect[cont]<min_value:
            min_value=vect[cont]
            pos=cont

    return pos

@jit(nopython=True, cache=True, fastmath=True)
def find_min_max(vect):
    min_value=min(vect)
    max_value=max(vect)
    return min_value,max_value

@jit(nopython=True, cache=True, fastmath=True)
def create_min_max_v_from_mesh(v0,v1,v2):
    N_points=v0.shape[0]
    meshXYZmin=np.zeros((N_points,3), dtype="double")
    meshXYZmax = np.zeros((N_points, 3), dtype="double")
    temp_V = np.zeros((3), dtype="double")
    for i in range(N_points):
        for j in range(3):
            temp_V[0] = v0[i, j]
            temp_V[1] = v1[i, j]
            temp_V[2] = v2[i, j]
            meshXYZmin[i,j],meshXYZmax[i,j] = find_min_max(temp_V)

    return meshXYZmin,meshXYZmax

@jit(nopython=True, cache=True, fastmath=True)
def find_pos_min_max_indices_conditioned(V_min,V_max,dimension,leq,geq):
    N=V_min.shape[0]
    ind=np.zeros((N), dtype="int64")
    pos=-1
    for i in range(N):
        if V_min[i,dimension]<=leq and V_max[i,dimension]>=geq:
            pos=pos+1
            ind[pos]=i

    if pos>=0:
        ind=ind[0:pos+1]
    else:
        ind = np.zeros((1), dtype="int64")
        ind[0] = -1

    return ind

@jit(nopython=True, cache=True, fastmath=True)
def find_pos_min_max_indices_conditioned_and_indicized(V_min,V_max,indices,dimension,leq,geq):
    N=len(indices)
    ind=np.zeros((N), dtype="int64")
    pos=-1
    for i in range(N):
        if (V_min[indices[i],dimension]<=leq) and (V_max[indices[i],dimension]>=geq):
            pos = pos + 1
            ind[pos]=i

    if pos >= 0:
        ind=ind[0:pos+1]
    else:
        ind = np.zeros((1), dtype="int64")
        ind[0] = -1

    return ind

@jit(nopython=True, cache=True, fastmath=True)
def find_ind_cross_list(v0,v1,v2,possibleCROSSLIST,gridCOxl,gridCOyl):
    N = len(possibleCROSSLIST)
    ind = np.zeros((N), dtype="int64")
    pos = -1
    for i in range(N):
        if (v0[possibleCROSSLIST[i],0]==gridCOxl and v1[possibleCROSSLIST[i],0]==gridCOyl) or \
            (v0[possibleCROSSLIST[i], 1] == gridCOxl and v1[possibleCROSSLIST[i], 1] == gridCOyl) or \
                (v0[possibleCROSSLIST[i], 2] == gridCOxl and v1[possibleCROSSLIST[i], 2] == gridCOyl):

            pos = pos + 1
            ind[pos] = i

    if pos >= 0:
        ind=ind[0:pos+1]
    else:
        ind = np.zeros((1), dtype="int64")
        ind[0] = -1
    return ind

@jit(nopython=True, cache=True, fastmath=True)
def find_first(V,val):
    N=len(V)
    pos=-1
    for i in range(N):
        if val==V[i]:
            pos=i
            break
    return pos

@jit(nopython=True, cache=True, fastmath=True)
def find_ind_to_keep_gridC0z(gridCOzCROSS,meshZmin,meshZmax):
    N=len(gridCOzCROSS)
    ind=np.zeros((N), dtype="int64")
    pos=-1
    for i in range(N):
        if (gridCOzCROSS[i]>=meshZmin-1e-12) and (gridCOzCROSS[i]<=meshZmax+1e-12):
            pos = pos + 1
            ind[pos]=i

    if pos >= 0:
        ind=ind[0:pos+1]
    else:
        ind = np.zeros((1), dtype="int64")
        ind[0] = -1

    return ind

@jit(nopython=True, cache=True, fastmath=True)
def find_ind_voxel_inside(V,val1,val2):
    N=len(V)
    ind=np.zeros((N), dtype="int64")
    pos=-1
    for i in range(N):
        if V[i]>val1 and V[i]<val2:
            pos=pos+1
            ind[pos]=i
    if pos >= 0:
        ind=ind[0:pos+1]
    else:
        ind = np.zeros((1), dtype="int64")
        ind[0] = -1

    return ind

@jit(nopython=True, cache=True, fastmath=True)
def find_all_equal_on_second_dimension(Mat2D,dimension,val):
    N=Mat2D.shape[0]
    ind = np.zeros((N), dtype="int64")
    pos = -1
    for i in range(N):
        if val==Mat2D[i,dimension]:
            pos = pos + 1
            ind[pos] = i
        if pos >= 0:
            ind = ind[0:pos + 1]
        else:
            ind = np.zeros((1), dtype="int64")
            ind[0] = -1

        return ind

@jit(nopython=True, cache=True, fastmath=True)
def CONVERT_meshformat(v0,v1,v2):

    vertices=np.vstack((v0, v1))
    vertices = np.vstack((vertices, v2))
    vertices = np.unique(vertices, axis=0)

    faces = np.zeros((v0.shape[0], 3), dtype="double")

    for loopF in range(v0.shape[0]):
        for loopV in range(3):

            if loopV==0:
                ind_1 = find_all_equal_on_second_dimension(vertices, 0, v0[loopF, 0])
                ind_2 = find_all_equal_on_second_dimension(vertices, 1, v0[loopF, 1])
                ind_3 = find_all_equal_on_second_dimension(vertices, 2, v0[loopF, 2])
            else:
                if loopV == 1:
                    ind_1 = find_all_equal_on_second_dimension(vertices, 0, v1[loopF, 0])
                    ind_2 = find_all_equal_on_second_dimension(vertices, 1, v1[loopF, 1])
                    ind_3 = find_all_equal_on_second_dimension(vertices, 2, v1[loopF, 2])
                else:
                    ind_1 = find_all_equal_on_second_dimension(vertices, 0, v2[loopF, 0])
                    ind_2 = find_all_equal_on_second_dimension(vertices, 1, v2[loopF, 1])
                    ind_3 = find_all_equal_on_second_dimension(vertices, 2, v2[loopF, 2])

            if ind_1[0]!=-1 and ind_2[0]!=-1 and ind_3[0]!=-1:
                vertref = np.intersect1d(ind_1, ind_2)
                vertref = np.intersect1d(vertref, ind_3)
                if len(vertref)>0:
                    faces[loopF, loopV] = vertref[0]

    return faces,vertices

@jit(nopython=True, cache=True, fastmath=True)
def COMPUTE_mesh_normals(v0,v1,v2):

    facetCOUNT = v0.shape[0]
    coordNORMALS = np.zeros((facetCOUNT,3), dtype="double")

    for loopFACE in range(facetCOUNT):

        # Find the coordinates for each vertex
        cornerA = v0[loopFACE, :]
        cornerB = v1[loopFACE, :]
        cornerC = v2[loopFACE, :]

        # Compute the vectors AB and AC
        AB = cornerB - cornerA
        AC = cornerC - cornerA

        # Determine the cross product AB x AC
        ABxAC = np.cross(AB, AC)

        # Normalise to give a unit vector
        ABxAC = ABxAC / LA.norm(ABxAC)
        coordNORMALS[loopFACE, :] = ABxAC

    return coordNORMALS

def voxel_intern(grid_x,grid_y,grid_z,v0_in,v1_in,v2_in,input_desc,case_perm):

    voxcountX = len(grid_x)
    voxcountY = len(grid_y)
    voxcountZ = len(grid_z)

    gridOUTPUT = np.full((voxcountX,voxcountY,voxcountZ), False, dtype=bool)

    Nnodes=v0_in.shape[0]
    assert case_perm in [0,1,2]
    if case_perm==0:
        meshXmin = input_desc['meshYmin']
        meshXmax = input_desc['meshYmax']
        meshYmin = input_desc['meshZmin']
        meshYmax = input_desc['meshZmax']
        meshZmin = input_desc['meshXmin']
        meshZmax = input_desc['meshXmax']
        v0 = np.zeros((Nnodes, 3), dtype="double")
        v1 = np.zeros((Nnodes, 3), dtype="double")
        v2 = np.zeros((Nnodes, 3), dtype="double")

        v0[:, 0]=  v0_in[:, 1]
        v1[:, 0] = v1_in[:, 1]
        v2[:, 0] = v2_in[:, 1]
        v0[:, 1] = v0_in[:, 2]
        v1[:, 1] = v1_in[:, 2]
        v2[:, 1] = v2_in[:, 2]
        v0[:, 2] = v0_in[:, 0]
        v1[:, 2] = v1_in[:, 0]
        v2[:, 2] = v2_in[:, 0]
    elif case_perm==1:
        meshXmin = input_desc['meshZmin']
        meshXmax = input_desc['meshZmax']
        meshYmin = input_desc['meshXmin']
        meshYmax = input_desc['meshXmax']
        meshZmin = input_desc['meshYmin']
        meshZmax = input_desc['meshYmax']
        v0 = np.zeros((Nnodes, 3), dtype="double")
        v1 = np.zeros((Nnodes, 3), dtype="double")
        v2 = np.zeros((Nnodes, 3), dtype="double")

        v0[:, 0] = v0_in[:, 2]
        v1[:, 0] = v1_in[:, 2]
        v2[:, 0] = v2_in[:, 2]
        v0[:, 1] = v0_in[:, 0]
        v1[:, 1] = v1_in[:, 0]
        v2[:, 1] = v2_in[:, 0]
        v0[:, 2] = v0_in[:, 1]
        v1[:, 2] = v1_in[:, 1]
        v2[:, 2] = v2_in[:, 1]
    else:
        assert case_perm==2
        meshXmin = input_desc['meshXmin']
        meshXmax = input_desc['meshXmax']
        meshYmin = input_desc['meshYmin']
        meshYmax = input_desc['meshYmax']
        meshZmin = input_desc['meshZmin']
        meshZmax = input_desc['meshZmax']
        v0 = v0_in
        v1 = v1_in
        v2 = v2_in

    # %Identify the min and max x,y coordinates (pixels) of the mesh:

    vect_temp=np.zeros((grid_x.shape[0]), dtype="double")
    for cont in range(grid_x.shape[0]):
        vect_temp[cont]=abs(grid_x[cont]-meshXmin)
    meshXminp = find_pos_min(vect_temp, grid_x.shape[0])

    for cont in range(grid_x.shape[0]):
        vect_temp[cont]=abs(grid_x[cont]-meshXmax)
    meshXmaxp = find_pos_min(vect_temp, grid_x.shape[0])

    vect_temp = np.zeros((grid_y.shape[0]), dtype="double")
    for cont in range(grid_y.shape[0]):
        vect_temp[cont] = abs(grid_y[cont] - meshYmin)

    meshYminp = find_pos_min(vect_temp, grid_y.shape[0])

    for cont in range(grid_y.shape[0]):
        vect_temp[cont] = abs(grid_y[cont] - meshYmax)

    meshYmaxp = find_pos_min(vect_temp, grid_y.shape[0])

    # %Make sure min < max for the mesh coordinates:
    if meshXminp > meshXmaxp:
        temp=meshXminp
        meshXminp=meshXmaxp
        meshXmaxp=temp
    if meshYminp > meshYmaxp:
        temp = meshYminp
        meshYminp=meshYmaxp
        meshYmaxp=temp

    # %Identify the min and max x,y,z coordinates of each facet:
    meshXYZmin,meshXYZmax = create_min_max_v_from_mesh(v0,v1,v2)

    # %======================================================
    # % VOXELISE THE MESH
    # %======================================================

    correctionLIST = np.zeros((0,2), dtype="int64")   #Prepare to record all rays that fail the voxelisation.  This array is built on-the-fly, but since
                                                               #it ought to be relatively small should not incur too much of a speed penalty.

    shift_div=4.37463724e-14 #to avoid division by 0

    #shift_div=0

    # % Loop through each x,y pixel.
    # % The mesh will be voxelised by passing rays in the z-direction through
    # % each x,y pixel, and finding the locations where the rays cross the mesh.
    for loopY in range (meshYminp,meshYmaxp+1):

        #   % - 1a - Find which mesh facets could possibly be crossed by the ray:
        ind_pcc = find_pos_min_max_indices_conditioned(meshXYZmin, meshXYZmax, 1, grid_y[loopY], grid_y[loopY] )

        if ind_pcc[0] != -1:
            possibleCROSSLISTy = ind_pcc

            for loopX in range (meshXminp,meshXmaxp+1):
                #     % - 1b - Find which mesh facets could possibly be crossed by the ray:
                ind_pos=find_pos_min_max_indices_conditioned_and_indicized(\
                        meshXYZmin, meshXYZmax,possibleCROSSLISTy ,0, grid_x[loopX], grid_x[loopX])

                if ind_pos[0]!=-1:

                    possibleCROSSLIST = possibleCROSSLISTy[ind_pos]

                    if len(possibleCROSSLIST)>0:
                        #Only continue the analysis if some nearby facets were actually identified
                        #  - 2 - For each facet, check if the ray really does cross the facet rather than just passing it close-by:
                        # GENERAL METHOD:
                        # A. Take each edge of the facet in turn.
                        # B. Find the position of the opposing vertex to that edge.
                        # C. Find the position of the ray relative to that edge.
                        # D. Check if ray is on the same side of the edge as the opposing vertex.
                        # E. If this is true for all three edges, then the ray definitely passes through the facet.
                        #
                        # NOTES:
                        # A. If a ray crosses exactly on a vertex:
                        #   a. If the surrounding facets have normal components pointing in the same (or opposite) direction as the ray then the face IS crossed.
                        #   b. Otherwise, add the ray to the correctionlist.

                        facetCROSSLIST = np.zeros((0), dtype="int64")  #Prepare to record all facets which are crossed by the ray.  This array is built on-the-fly, but since
                        # it ought to be relatively small (typically a list of <10) should not incur too much of a speed penalty.

                        #----------
                        # - 1 - Check for crossed vertices:
                        #---------

                        #Find which mesh facets contain a vertex which is crossed by the ray:
                        ind_ccl=find_ind_cross_list(v0,v1,v2, possibleCROSSLIST, grid_x[loopX], grid_y[loopY])
                        if ind_ccl[0] != -1:
                            vertexCROSSLIST=possibleCROSSLIST[ind_ccl]
                            checkindex = np.zeros((len(vertexCROSSLIST)), dtype="int64")

                            continue_cycle=True
                            while continue_cycle==True:
                                vertexindex=find_first(checkindex, 0)
                                if vertexindex==-1:
                                    continue_cycle=False
                                else:
                                    vertexindex = find_first(checkindex, -1)
                                    temp_faces, temp_vertices = CONVERT_meshformat(v0[vertexCROSSLIST,:],\
                                                                v1[vertexCROSSLIST,:],v2[vertexCROSSLIST,:])
                                    adjacentindex = np.isin(temp_faces, temp_faces[vertexindex,:])
                                    coN=np.zeros((adjacentindex.shape[0],3), dtype="double")
                                    p_in=-1
                                    for caux in range(adjacentindex.shape[0]):
                                        if any(adjacentindex[caux,:])==True:
                                            checkindex[caux]=1
                                            p_in=p_in+1
                                            coN[p_in,:]= COMPUTE_mesh_normals(v0[vertexCROSSLIST[caux],:,:], \
                                                v1[vertexCROSSLIST[caux],:,:],v2[vertexCROSSLIST[caux], :, :] )

                                    if p_in!=-1:
                                        minco,maxc0=find_min_max(coN[:,2], coN.shape[0])
                                        if maxc0<0 or minco>0:
                                            facetCROSSLIST    = np.hstack((facetCROSSLIST, vertexCROSSLIST[vertexindex]))
                                    else:
                                        possibleCROSSLIST = np.zeros((0), dtype="int64")
                                        correctionLIST    = np.vstack((correctionLIST, np.hstack((loopX, loopY))))
                                        checkindex = np.ones((len(vertexCROSSLIST)), dtype="int64")

                    #----------
                    # - 2 - Check for crossed facets:
                    #----------
                    Npc=len(possibleCROSSLIST)
                    if Npc>0:  #Only continue the analysis if some nearby facets were actually identified
                        for i in range(Npc):
                            loopCHECKFACET= possibleCROSSLIST[i]

                            if loopCHECKFACET!=-1:

                                #Check if ray crosses the facet.  This method is much (>>10 times) faster than using the built-in function 'inpolygon'.
                                #Taking each edge of the facet in turn, check if the ray is on the same side as the opposing vertex.

                                Y1predicted = v1[loopCHECKFACET,1] - ((v1[loopCHECKFACET,1] -v2[loopCHECKFACET,1] ) * (v1[loopCHECKFACET,0] -v0[loopCHECKFACET,0] )/(shift_div+v1[loopCHECKFACET,0] -v2[loopCHECKFACET,0] ))
                                YRpredicted = v1[loopCHECKFACET,1]  - ((v1[loopCHECKFACET,1] -v2[loopCHECKFACET,1] ) * (v1[loopCHECKFACET,0] -grid_x[loopX])/(shift_div+v1[loopCHECKFACET,0]-v2[loopCHECKFACET,0]))

                                if (Y1predicted > v0[loopCHECKFACET,1] and YRpredicted > grid_y[loopY]) or (Y1predicted < v0[loopCHECKFACET,1] and YRpredicted < grid_y[loopY]):
                                    #The ray is on the same side of the 2-3 edge as the 1st vertex.

                                    Y2predicted = v2[loopCHECKFACET,1] - ((v2[loopCHECKFACET,1]-v0[loopCHECKFACET,1]) * (v2[loopCHECKFACET,0]-v1[loopCHECKFACET,0])/(shift_div+v2[loopCHECKFACET,0]-v0[loopCHECKFACET,0]))
                                    YRpredicted = v2[loopCHECKFACET,1] - ((v2[loopCHECKFACET,1]-v0[loopCHECKFACET,1]) * (v2[loopCHECKFACET,0]-grid_x[loopX])/(shift_div+v2[loopCHECKFACET,0]-v0[loopCHECKFACET,0]))

                                    if (Y2predicted > v1[loopCHECKFACET,1] and YRpredicted > grid_y[loopY]) or (Y2predicted < v1[loopCHECKFACET,1] and YRpredicted < grid_y[loopY]):
                                        #The ray is on the same side of the 3-1 edge as the 2nd vertex.
                                        Y3predicted = v0[loopCHECKFACET,1] - ((v0[loopCHECKFACET,1]-v1[loopCHECKFACET,1]) * (v0[loopCHECKFACET,0]-v2[loopCHECKFACET,0])/(shift_div+v0[loopCHECKFACET,0]-v1[loopCHECKFACET,0]))
                                        YRpredicted = v0[loopCHECKFACET,1] - ((v0[loopCHECKFACET,1]-v1[loopCHECKFACET,1]) * (v0[loopCHECKFACET,0]-grid_x[loopX])/(shift_div+v0[loopCHECKFACET,0]-v1[loopCHECKFACET,0]))

                                        if (Y3predicted > v2[loopCHECKFACET,1] and YRpredicted > grid_y[loopY]) or (Y3predicted < v2[loopCHECKFACET,1] and YRpredicted < grid_y[loopY]):
                                            # The ray is on the same side of the 1-2 edge as the 3rd vertex.
                                            #The ray passes through the facet since it is on the correct side of all 3 edges
                                            facetCROSSLIST = np.hstack((facetCROSSLIST, loopCHECKFACET))

                    #----------
                    # - 3 - Find the z coordinate of the locations where the ray crosses each facet or vertex:
                    #----------

                    Nfc = len(facetCROSSLIST)
                    grid_zCROSS = np.zeros((Nfc), dtype="double")
                    for i in range(Nfc):
                        loopFINDZ=facetCROSSLIST[i]

                        #  METHOD:
                        # 1. Define the equation describing the plane of the facet.  For a
                        # more detailed outline of the maths, see:
                        # http://local.wasp.uwa.edu.au/~pbourke/geometry/planeeq/
                        #    Ax + By + Cz + D = 0
                        #    where  A = y1 (z2 - z3) + y2 (z3 - z1) + y3 (z1 - z2)
                        #           B = z1 (x2 - x3) + z2 (x3 - x1) + z3 (x1 - x2)
                        #           C = x1 (y2 - y3) + x2 (y3 - y1) + x3 (y1 - y2)
                        #           D = - x1 (y2 z3 - y3 z2) - x2 (y3 z1 - y1 z3) - x3 (y1 z2 - y2 z1)
                        # 2. For the x and y coordinates of the ray, solve these equations to find the z coordinate in this plane.

                        planecoA = v0[loopFINDZ,1]*(v1[loopFINDZ,2]-v2[loopFINDZ,2]) + v1[loopFINDZ,1]*(v2[loopFINDZ,2]-v0[loopFINDZ,2]) + v2[loopFINDZ,1]*(v0[loopFINDZ,2]-v1[loopFINDZ,2])
                        planecoB = v0[loopFINDZ,2]*(v1[loopFINDZ,0]-v2[loopFINDZ,0]) + v1[loopFINDZ,2]*(v2[loopFINDZ,0]-v0[loopFINDZ,0]) + v2[loopFINDZ,2]*(v0[loopFINDZ,0]-v1[loopFINDZ,0])
                        planecoC = v0[loopFINDZ,0]*(v1[loopFINDZ,1]-v2[loopFINDZ,1]) + v1[loopFINDZ,0]*(v2[loopFINDZ,1]-v0[loopFINDZ,1]) + v2[loopFINDZ,0]*(v0[loopFINDZ,1]-v1[loopFINDZ,1])
                        planecoD = - v0[loopFINDZ,0]*(v1[loopFINDZ,1]*v2[loopFINDZ,2]-v2[loopFINDZ,1]*v1[loopFINDZ,2]) - v1[loopFINDZ,0]*(v2[loopFINDZ,1]*v0[loopFINDZ,2]-v0[loopFINDZ,1]*v2[loopFINDZ,2]) \
                                   - v2[loopFINDZ,0]*(v0[loopFINDZ,1]*v1[loopFINDZ,2]-v1[loopFINDZ,1]*v0[loopFINDZ,2])

                        if abs(planecoC) < 1e-14:
                            planecoC=0.0

                        grid_zCROSS[i] = (- planecoD - planecoA*grid_x[loopX] - planecoB*grid_y[loopY]) / (shift_div+planecoC)

                    # Remove values of grid_zCROSS which are outside of the mesh limits (including a 1e-12 margin for error).
                    ind_pos_keep=find_ind_to_keep_gridC0z(grid_zCROSS, meshZmin, meshZmax)
                    if ind_pos_keep[0] != -1:
                        grid_zCROSS=grid_zCROSS[ind_pos_keep]

                        #Round grid_zCROSS to remove any rounding errors, and take only the unique values:
                        grid_zCROSS = np.round(grid_zCROSS*1e12)/1e12
                        grid_zCROSS = np.unique(grid_zCROSS)

                        #----------
                        # - 4 - Label as being inside the mesh all the voxels that the ray passes through after crossing one facet before crossing another facet:
                        #----------
                        if np.remainder(len(grid_zCROSS), 2)==0: #Only rays which cross an even number of facets are voxelised

                            for loopASSIGN in range(np.int64(np.floor(len(grid_zCROSS)/2))):
                                voxelsINSIDE = find_ind_voxel_inside(grid_z,grid_zCROSS[2*(loopASSIGN+1)-2],grid_zCROSS[2*(loopASSIGN+1)-1])
                                if voxelsINSIDE[0] != -1:
                                    gridOUTPUT[loopX,loopY,voxelsINSIDE] = True

                        else:
                            if len(grid_zCROSS)>0:  # Remaining rays which meet the mesh in some way are not voxelised, but are labelled for correction later.
                                correctionLIST=np.vstack((correctionLIST, np.hstack((loopX, loopY))))

        # ======================================================
        #  USE INTERPOLATION TO FILL IN THE RAYS WHICH COULD NOT BE VOXELISED
        # ======================================================
        # For rays where the voxelisation did not give a clear result, the ray is
        # computed by interpolating from the surrounding rays.

    countCORRECTIONLIST = correctionLIST.shape[0]

    if countCORRECTIONLIST>0:
    
        # If necessary, add a one-pixel border around the x and y edges of the
        # array.  This prevents an error if the code tries to interpolate a ray at
        # the edge of the x,y grid.
        if min(correctionLIST[:,0])==0 or max(correctionLIST[:,0])==(len(grid_x)-1) or \
                min(correctionLIST[:,1])==0 or max(correctionLIST[:,1])==(len(grid_y)-1):

            temp = np.hstack((np.full((voxcountX, 1, voxcountZ), False, dtype=bool),gridOUTPUT))
            temp = np.hstack((temp, np.full((voxcountX, 1, voxcountZ), False, dtype=bool)))
            temp = np.vstack((np.full((1, voxcountY+2, voxcountZ), False, dtype=bool),temp))
            gridOUTPUT = np.vstack((temp, np.full((1, voxcountY+2, voxcountZ), False, dtype=bool)))
            correctionLIST = correctionLIST + np.ones((countCORRECTIONLIST,2), dtype="int64")

  
        for loopC in range(countCORRECTIONLIST):
            for cz in range(gridOUTPUT.shape[2]):
                voxelsforcorrection = 0
                if gridOUTPUT[correctionLIST[loopC,0]-1,correctionLIST[loopC,1]-1,cz]==True:
                    voxelsforcorrection=voxelsforcorrection+1
                if gridOUTPUT[correctionLIST[loopC,0]-1,correctionLIST[loopC,1],cz]==True:
                    voxelsforcorrection=voxelsforcorrection+1
                if gridOUTPUT[correctionLIST[loopC, 0] - 1, correctionLIST[loopC, 1]+1, cz] == True:
                    voxelsforcorrection=voxelsforcorrection+1
                if gridOUTPUT[correctionLIST[loopC, 0] , correctionLIST[loopC, 1] , cz] == True:
                    voxelsforcorrection=voxelsforcorrection+1
                if gridOUTPUT[correctionLIST[loopC, 0] , correctionLIST[loopC, 1] + 1, cz] == True:
                    voxelsforcorrection=voxelsforcorrection+1
                if gridOUTPUT[correctionLIST[loopC, 0] + 1, correctionLIST[loopC, 1] - 1, cz] == True:
                    voxelsforcorrection=voxelsforcorrection+1
                if gridOUTPUT[correctionLIST[loopC, 0] + 1, correctionLIST[loopC, 1] , cz] == True:
                    voxelsforcorrection=voxelsforcorrection+1
                if gridOUTPUT[correctionLIST[loopC, 0] + 1, correctionLIST[loopC, 1] + 1, cz] == True:
                    voxelsforcorrection=voxelsforcorrection+1

                if voxelsforcorrection>3:
                    gridOUTPUT[correctionLIST[loopC, 0], correctionLIST[loopC, 1], voxelsforcorrection] = True

        #Remove the one-pixel border surrounding the array, if this was added previously
        Ntx = gridOUTPUT.shape[0]
        Nty = gridOUTPUT.shape[1]

        if Ntx>len(grid_x) or Nty>len(grid_y):
            gridOUTPUT = gridOUTPUT[1:Ntx-1,1:Nty-1,:]

    return gridOUTPUT
# %==========================================================================