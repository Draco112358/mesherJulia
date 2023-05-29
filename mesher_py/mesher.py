try:
  import unzip_requirements
except ImportError:
  pass
# handler.py
from stl.mesh import Mesh
# import vtkplotlib as vpl
import os
# import csv
import sys
import math
import numpy as np
# import scipy
from .voxelizator import voxelize, solve_overlapping
import json


def find_mins_maxs(mesh_object: Mesh):
    assert type(mesh_object) == Mesh
    minx = mesh_object.x.min()
    maxx = mesh_object.x.max()
    miny = mesh_object.y.min()
    maxy = mesh_object.y.max()
    minz = mesh_object.z.min()
    maxz = mesh_object.z.max()
    return minx, maxx, miny, maxy, minz, maxz


def find_box_dimensions(dict_meshes: dict):
    global_min_x, global_min_y, global_min_z = sys.maxsize, sys.maxsize, sys.maxsize
    global_max_x, global_max_y, global_max_z = -sys.maxsize, -sys.maxsize, -sys.maxsize

    for key in dict_meshes:
        ms = dict_meshes[key]
        minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(ms)
        global_min_x = min(global_min_x, minx)
        global_min_y = min(global_min_y, miny)
        global_min_z = min(global_min_z, minz)
        global_max_x = max(global_max_x, maxx)
        global_max_y = max(global_max_y, maxy)
        global_max_z = max(global_max_z, maxz)

    keeper_object = {}
    keeper_object['meshXmin'] = global_min_x
    keeper_object['meshXmax'] = global_max_x
    keeper_object['meshYmin'] = global_min_y
    keeper_object['meshYmax'] = global_max_y
    keeper_object['meshZmin'] = global_min_z
    keeper_object['meshZmax'] = global_max_z

    w = keeper_object['meshXmax'] - keeper_object['meshXmin']
    l = keeper_object['meshYmax'] - keeper_object['meshYmin']
    h = keeper_object['meshZmax'] - keeper_object['meshZmin']

    return w, l, h, keeper_object


def find_sizes(number_of_cells_x: int, number_of_cells_y: int, number_of_cells_z: int, geometry_descriptor: dict):
    
    assert type(number_of_cells_x) == int
    assert type(number_of_cells_y) == int
    assert type(number_of_cells_z) == int
    assert type(geometry_descriptor) == dict
    assert len(geometry_descriptor) == 6

    # minimum_vertex_coordinates = [geometry_descriptor['meshXmin'] * 1e-3, geometry_descriptor['meshYmin'] * 1e-3,
    #           geometry_descriptor['meshZmin'] * 1e-3]
    # # max_v = [minmax.meshXmax minmax.meshYmax minmax.meshZmax]*1e-3;
    xv = np.linspace(geometry_descriptor['meshXmin'] * 1e-3, geometry_descriptor['meshXmax'] * 1e-3,
                     number_of_cells_x + 1)
    yv = np.linspace(geometry_descriptor['meshYmin'] * 1e-3, geometry_descriptor['meshYmax'] * 1e-3,
                     number_of_cells_y + 1)
    zv = np.linspace(geometry_descriptor['meshZmin'] * 1e-3, geometry_descriptor['meshZmax'] * 1e-3,
                     number_of_cells_z + 1)

    return abs(xv[2] - xv[1]), abs(yv[2] - yv[1]), abs(zv[2] - zv[1])#, minimum_vertex_coordinates

def dump_json_data(filename,n_materials,o_x,o_y,o_z,cs_x,cs_y,cs_z,nc_x,nc_y,nc_z,matr,id_to_material):
    
    #print("Serialization to:",filename)
    assert(isinstance(matr, np.ndarray))
    assert type(cs_x) == np.float64
    assert type(cs_y) == np.float64
    assert type(cs_z) == np.float64
    assert type(o_x) == np.float64
    assert type(o_y) == np.float64
    assert type(o_z) == np.float64

    origin = {"origin_x":float(o_x),"origin_y":float(o_y),"origin_z":float(o_z)}
    
    n_cells = {"n_cells_x":float(nc_x),"n_cells_y":float(nc_y),"n_cells_z":float(nc_z)}

    cell_size = {"cell_size_x":float(cs_x),"cell_size_y":float(cs_y),"cell_size_z":float(cs_z)}

    materials = {}
    for element in id_to_material:
        materials[element]=id_to_material[element]
        
    mesher_matrices_dict = {}
    
    count = 0
    
    for matrix in matr.tolist():
        assert count in id_to_material 
        mesher_matrices_dict[id_to_material[count]] = matrix
        count += 1
    assert count == n_materials
    json_dict = {"n_materials" :n_materials,"materials": materials, "origin": origin , "cell_size": cell_size, "n_cells" : n_cells, "mesher_matrices": mesher_matrices_dict}
    return {
        'statusCode': 200,
        'body': json_dict
        }

        
def doMeshing(event, context):

    meshes = {}
    print(event)
    
    dictData = json.loads(event['body'])
    type(dictData['STLList'])==list
    for geometry in dictData['STLList']:
        assert type(geometry)==dict
        mesh_id = geometry['material']
        mesh_stl = geometry['STL']
        assert mesh_id not in meshes
        with open("/tmp/stl.temp", "w") as write_file:
            write_file.write(mesh_stl)
        mesh_stl_converted = Mesh.from_file("/tmp/stl.temp")
        assert type(mesh_stl_converted)==Mesh
        meshes[mesh_id] = mesh_stl_converted
        os.remove("/tmp/stl.temp") 

    
    geometry_x_bound, geometry_y_bound, geometry_z_bound, geometry_data_object = find_box_dimensions(dict_meshes=meshes)
    

    # grids grain
    assert type(dictData['quantum'])==list
    quantum_x, quantum_y, quantum_z = dictData['quantum']

    # quantum_x, quantum_y, quantum_z = 1, 1e-2, 1e-2 #per Test 1
    # # quantum_x, quantum_y, quantum_z = 1e-1, 1, 1e-2  # per Test 2
    # # quantum_x, quantum_y, quantum_z = 1e-1, 1e-1, 1e-2  # per Test 3
    # # quantum_x, quantum_y, quantum_z = 2, 1, 1e-2  # per Test 4
    # # quantum_x, quantum_y, quantum_z = 1, 1, 1e-2  # per Test 5

    #print("QUANTA:",quantum_x, quantum_y, quantum_z)

    n_of_cells_x = math.ceil(geometry_x_bound / quantum_x)
    n_of_cells_y = math.ceil(geometry_y_bound / quantum_y)
    n_of_cells_z = math.ceil(geometry_z_bound / quantum_z)
    
    #print("GRID:",n_of_cells_x, n_of_cells_y, n_of_cells_z)
    
    cell_size_x, cell_size_y, cell_size_z = find_sizes(number_of_cells_x=n_of_cells_x,
                                                              number_of_cells_y=n_of_cells_y,
                                                              number_of_cells_z=n_of_cells_z,
                                                              geometry_descriptor=geometry_data_object)
    
    precision = 0.1
    #print("CELL SIZE AFTER ADJUSTEMENTS:",(cell_size_x), (cell_size_y), (cell_size_z))
    # if __debug__:
        
    #     for size,quantum in zip([cell_size_x,cell_size_y,cell_size_z],[quantum_x,quantum_y,quantum_z]):
    #         print(abs(size*(1/precision) - quantum),precision)
    #         assert abs(size*(1/precision) - quantum)<=precision
            

    
    n_materials = len(dictData['STLList'])
    
    mesher_output = np.full((n_materials,n_of_cells_x, n_of_cells_y, n_of_cells_z), False, dtype=bool)

    mapping_ids_to_materials = {}

    counter_stl_files = 0
    for mesh_id in meshes:
        
        assert type(meshes[mesh_id])==Mesh
        #print("voxeling",mesh_id)

        mesher_output[counter_stl_files,:,:,:] = voxelize(cells_on_x=n_of_cells_x, cells_on_y=n_of_cells_y, cells_on_z=n_of_cells_z,
                                      meshXYZ=meshes[mesh_id], geometry_desc=geometry_data_object)

        mapping_ids_to_materials[counter_stl_files]=mesh_id
        counter_stl_files+=1

    id_mats_keep=np.zeros((n_materials), dtype="int64")
    
    id_mats_keep[0]=0
    
    solve_overlapping(n_cells_x=n_of_cells_x, n_cells_y=n_of_cells_y, n_cells_z=n_of_cells_z, num_materials=n_materials, id_mat_keep=id_mats_keep, output_meshing=mesher_output)
        



    origin_x = geometry_data_object['meshXmin']*1e-3
    origin_y = geometry_data_object['meshYmin']*1e-3
    origin_z = geometry_data_object['meshZmin']*1e-3


    assert(isinstance(mesher_output, np.ndarray))
    assert type(cell_size_x) == np.float64
    assert type(cell_size_y) == np.float64
    assert type(cell_size_z) == np.float64
    assert type(origin_x) == np.float64
    assert type(origin_y) == np.float64
    assert type(origin_z) == np.float64


    # Writing to data.json
    json_file = "outputMesher.json"
    
    return dump_json_data(filename=json_file,n_materials=counter_stl_files,o_x=origin_x,o_y=origin_y,o_z=origin_z,\
            cs_x=cell_size_x,cs_y=cell_size_y,cs_z=cell_size_z,\
            nc_x=n_of_cells_x,nc_y=n_of_cells_y,nc_z=n_of_cells_z,\
            matr=mesher_output,id_to_material=mapping_ids_to_materials) 
    


