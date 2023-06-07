
function solve_overlapping(n_cells_x,n_cells_y,n_cells_z,num_materials,id_mat_keep,output_meshing)

    for c1 in range(1, n_cells_x)
        for c2 in range(1, n_cells_y)
            for c3 in range(1, n_cells_z)
                for k in range(1, num_materials)
                    if output_meshing[k,c1, c2, c3] == true
                        for k2 in range(num_materials)
                            if output_meshing[k2,c1, c2, c3] == true && k!=k2
                                if k in id_mat_keep
                                    output_meshing[k2, c1, c2, c3] = false
                                else
                                    output_meshing[k, c1, c2, c3] = false
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function merge_the_3_grids(voxcountX, voxcountY, voxcountZ,gridOUTPUT1,gridOUTPUT2,gridOUTPUT3,gridOUTPUT)
    for c1 in range(1, voxcountX)
        for c2 in range(1, voxcountY)
            for c3 in range(1,voxcountZ)
                cont_true = 0
                if gridOUTPUT1[c2, c3, c1] == true
                    cont_true = cont_true + 1
                end
                if gridOUTPUT2[c3, c1, c2] == true
                    cont_true = cont_true + 1
                end
                if gridOUTPUT3[c1, c2, c3] == true
                    cont_true = cont_true + 1
                end
                if cont_true > 1
                    gridOUTPUT[c1, c2, c3] = true
                end
            end
        end
    end
end

function voxelize(cells_on_x,cells_on_y,cells_on_z,meshXYZ, geometry_desc)
       
    meshXmin = geometry_desc["meshXmin"]
    meshXmax = geometry_desc["meshXmax"]
    meshYmin = geometry_desc["meshYmin"]
    meshYmax = geometry_desc["meshYmax"]
    meshZmin = geometry_desc["meshZmin"]
    meshZmax = geometry_desc["meshZmax"]


    if cells_on_x==1
        #If gridX is a single integer (rather than a vector) and is equal to 1
        gridCOx   = (meshXmin+meshXmax)/2
    else
        #If gridX is a single integer (rather than a vector) then automatically create the list of x coordinates
        voxwidth  = (meshXmax-meshXmin)/(cells_on_x+1/2)
        gridCOx   = range(meshXmin+voxwidth/2, meshXmax-voxwidth/2, step=voxwidth)
    end
    
    if cells_on_y==12
        #If gridX is a single integer (rather than a vector) and is equal to 1
        gridCOy   = (meshYmin+meshYmax)/2
    else
        #If gridX is a single integer (rather than a vector) then automatically create the list of x coordinates
        voxwidth  = (meshYmax-meshYmin)/(cells_on_y+1/2)
        gridCOy   = range(meshYmin+voxwidth/2, meshYmax-voxwidth/2, step=voxwidth)
    end
    
    if cells_on_z==1
        #If gridX is a single integer (rather than a vector) and is equal to 1
        gridCOz   = (meshZmin+meshZmax)/2
    else
        #If gridX is a single integer (rather than a vector) then automatically create the list of x coordinates
        voxwidth  = (meshZmax-meshZmin)/(cells_on_z+1/2)
        gridCOz   = range(meshZmin+voxwidth/2, meshZmax-voxwidth/2, step=voxwidth)
    end

    # Check that the output grid is large enough to cover the mesh:
    var_check=0
    if (min(gridCOx)>meshXmin || max(gridCOx)<meshXmax)
        var_check = 1
        gridcheckX = 0
        if (min(gridCOx)>meshXmin)
            gridCOx=hcat((meshXmin, gridCOx))
            gridcheckX = gridcheckX+1
        end
        if (max(gridCOx)<meshXmax)
            gridCOx = hcat((gridCOx,meshXmax))
            gridcheckX = gridcheckX+2
        end
    else
        if (min(gridCOy)>meshYmin || max(gridCOy)<meshYmax)
            var_check = 2
            gridcheckY = 0
            if min(gridCOy)>meshYmin
                gridCOy = hcat((meshYmin, gridCOy))
                gridcheckY = gridcheckY+1
            end
            if max(gridCOy)<meshYmax
                gridCOy = hcat((gridCOy, meshYmax))
                gridcheckY = gridcheckY+2
            end
        else
            if (min(gridCOz)>meshZmin || max(gridCOz)<meshZmax)
                var_check = 3
                gridcheckZ = 0
                if (min(gridCOz)>meshZmin)
                    gridCOz = hcat((meshZmin, gridCOz))
                    gridcheckZ = gridcheckZ+1
                end
                if (max(gridCOz)<meshZmax)
                    gridCOz = hcat((gridCOz, meshZmax))
                    gridcheckZ = gridcheckZ+2
                end
            end
        end
    end

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

    gridOUTPUT = fill(false, (voxcountX, voxcountY, voxcountZ))
    merge_the_3_grids(voxcountX, voxcountY, voxcountZ, gridOUTPUT1, gridOUTPUT2, gridOUTPUT3,gridOUTPUT)

    # %======================================================
    # % RETURN THE OUTPUT GRID TO THE SIZE REQUIRED BY THE USER (IF IT WAS CHANGED EARLIER)
    # %======================================================
    # match var_check:
    if var_check == 1
        if gridcheckX == 1
            gridOUTPUT=gridOUTPUT[2:voxcountX+1,:,:]
            gridCOx    = gridCOx[2:voxcountX+1]
        elseif gridcheckX == 2
            gridOUTPUT = gridOUTPUT[1:voxcountX, :, :]
            gridCOx = gridCOx[1:voxcountX]
        else
            gridOUTPUT = gridOUTPUT[2:voxcountX, :, :]
            gridCOx = gridCOx[2:voxcountX]
        end
    elseif var_check == 2
            if gridcheckY == 1
                gridOUTPUT = gridOUTPUT[:,2:voxcountY+1, :]
                gridCOy = gridCOy[2:voxcountY+1]
            elseif gridcheckY == 2
                gridOUTPUT = gridOUTPUT[:, 1:voxcountY, :]
                gridCOy = gridCOy[1:voxcountY]
            else
                gridOUTPUT = gridOUTPUT[:, 2:voxcountY, :]
                gridCOy = gridCOy[2:voxcountY]
            end
    else
        if gridcheckZ == 1
            gridOUTPUT = gridOUTPUT[:, :, 2:voxcountZ+1]
            gridCOz = gridCOz[2:voxcountZ+1]
        elseif gridcheckZ == 2
            gridOUTPUT = gridOUTPUT[:, :, 1:voxcountZ]
            gridCOz = gridCOz[1:voxcountZ]
        else
            gridOUTPUT = gridOUTPUT[:, :, 2:voxcountZ]
            gridCOz = gridCOz[2:voxcountZ]
        end
    end

    #if needed also gridCOx, gridCOy and gridCOz can be given in output
    return gridOUTPUT
end
