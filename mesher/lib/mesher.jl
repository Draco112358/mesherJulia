include("voxelizator.jl")
include("voxelize_internal.jl")
using MeshIO
using FileIO
using Meshes
using MeshBridge

function find_mins_maxs(mesh_object::Mesh)
    bb = boundingbox(mesh_object)
    #@assert mesh_object isa Mesh
    minx =  coordinates(minimum(bb))[1]
    maxx = coordinates(maximum(bb))[1]
    miny = coordinates(minimum(bb))[2]
    maxy = coordinates(maximum(bb))[2]
    minz = coordinates(minimum(bb))[3]
    maxz = coordinates(maximum(bb))[3]
    return minx, maxx, miny, maxy, minz, maxz
end


function find_box_dimensions(dict_meshes::Dict)
    global_min_x, global_min_y, global_min_z = prevfloat(typemax(Float64)), prevfloat(typemax(Float64)), prevfloat(typemax(Float64))
    global_max_x, global_max_y, global_max_z = -prevfloat(typemax(Float64)), -prevfloat(typemax(Float64)), -prevfloat(typemax(Float64))

    for (key,value) in dict_meshes
        #println(dict_meshes)
        minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(value)
        global_min_x = min(global_min_x, minx)
        global_min_y = min(global_min_y, miny)
        global_min_z = min(global_min_z, minz)
        global_max_x = max(global_max_x, maxx)
        global_max_y = max(global_max_y, maxy)
        global_max_z = max(global_max_z, maxz)
    end

    keeper_object = Dict()
    keeper_object["meshXmin"] = global_min_x
    keeper_object["meshXmax"] = global_max_x
    keeper_object["meshYmin"] = global_min_y
    keeper_object["meshYmax"] = global_max_y
    keeper_object["meshZmin"] = global_min_z
    keeper_object["meshZmax"] = global_max_z

    w = keeper_object["meshXmax"] - keeper_object["meshXmin"]
    l = keeper_object["meshYmax"] - keeper_object["meshYmin"]
    h = keeper_object["meshZmax"] - keeper_object["meshZmin"]

    return w, l, h, keeper_object
end


function find_sizes(number_of_cells_x::Int, number_of_cells_y::Int, number_of_cells_z::Int, geometry_descriptor::Dict)
    
    @assert number_of_cells_x isa Int
    @assert number_of_cells_y isa Int
    @assert number_of_cells_z isa Int
    @assert geometry_descriptor isa Dict
    @assert length(geometry_descriptor)==6

    # minimum_vertex_coordinates = [geometry_descriptor['meshXmin'] * 1e-3, geometry_descriptor['meshYmin'] * 1e-3,
    #           geometry_descriptor['meshZmin'] * 1e-3]
    # # max_v = [minmax.meshXmax minmax.meshYmax minmax.meshZmax]*1e-3;
    xv = LinRange(geometry_descriptor["meshXmin"] * 1e-3, geometry_descriptor["meshXmax"] * 1e-3,
                     number_of_cells_x + 1)
    yv = LinRange(geometry_descriptor["meshYmin"] * 1e-3, geometry_descriptor["meshYmax"] * 1e-3,
                     number_of_cells_y + 1)
    zv = LinRange(geometry_descriptor["meshZmin"] * 1e-3, geometry_descriptor["meshZmax"] * 1e-3,
                     number_of_cells_z + 1)

    return abs(xv[3] - xv[2]), abs(yv[3] - yv[2]), abs(zv[3] - zv[2])#, minimum_vertex_coordinates
end

function dump_json_data(filename,n_materials,o_x::Float64,o_y::Float64,o_z::Float64,cs_x::Float64,cs_y::Float64,cs_z::Float64,nc_x,nc_y,nc_z,matr,id_to_material)
    
    #print("Serialization to:",filename)
    @assert cs_x isa Float64
    @assert cs_y isa Float64
    @assert cs_z isa Float64
    @assert o_x isa Float64
    @assert o_y isa Float64
    @assert o_z isa Float64

    origin = Dict("origin_x" => o_x, "origin_y" => o_y, "origin_z" => o_z)
    
    n_cells = Dict("n_cells_x" => parse(Float64, nc_x),"n_cells_y"=> parse(Float64, nc_y),"n_cells_z"=> parse(Float64, nc_z))

    cell_size = Dict("cell_size_x" => cs_x,"cell_size_y" => cs_y,"cell_size_z" => cs_z)

    materials = Dict()
    for element in id_to_material
        materials[element]=id_to_material[element]
    end
        
    mesher_matrices_dict = Dict()
    
    count = 1
    
    for matrix in matr.tolist()
        @assert count in id_to_material 
        mesher_matrices_dict[id_to_material[count]] = matrix
        count += 1
    end
    @assert count == n_materials+1
    json_dict = Dict("n_materials" => n_materials,"materials" => materials, "origin" => origin , "cell_size" => cell_size, "n_cells" => n_cells, "mesher_matrices" => mesher_matrices_dict)
    return json_dict
end

        
function doMeshing(dictData::Dict)

    meshes = Dict()
    for geometry in Array{Any}(dictData["STLList"])
        #@assert geometry isa Dict
        mesh_id = geometry["material"]
        mesh_stl = geometry["STL"]
        #@assert mesh_id not in meshes
        open("/tmp/stl.stl", "w") do write_file
            write(write_file,mesh_stl)
        end
        mesh_stl = load("/tmp/stl.stl")
        mesh_stl_converted = convert(Meshes.Mesh, mesh_stl)
        bb = boundingbox(mesh_stl_converted)
        println(coordinates(minimum(bb)))
        println(coordinates(maximum(bb)))
    
        #mesh_stl_converted = Meshes.Polytope(3,3,mesh_stl)
        #@assert mesh_stl_converted isa Mesh
        meshes[mesh_id] = mesh_stl_converted
        Base.Filesystem.rm("/tmp/stl.stl", force=true)
    end

    
    geometry_x_bound, geometry_y_bound, geometry_z_bound, geometry_data_object = find_box_dimensions(meshes)
    

    # grids grain
    # assert type(dictData['quantum'])==list
    quantum_x, quantum_y, quantum_z = dictData["quantum"]

    # quantum_x, quantum_y, quantum_z = 1, 1e-2, 1e-2 #per Test 1
    # # quantum_x, quantum_y, quantum_z = 1e-1, 1, 1e-2  # per Test 2
    # # quantum_x, quantum_y, quantum_z = 1e-1, 1e-1, 1e-2  # per Test 3
    # # quantum_x, quantum_y, quantum_z = 2, 1, 1e-2  # per Test 4
    # # quantum_x, quantum_y, quantum_z = 1, 1, 1e-2  # per Test 5

    #print("QUANTA:",quantum_x, quantum_y, quantum_z)

    n_of_cells_x = ceil(Int64, geometry_x_bound / quantum_x)
    n_of_cells_y = ceil(Int64, geometry_y_bound / quantum_y)
    n_of_cells_z = ceil(Int64, geometry_z_bound / quantum_z)
    
    #print("GRID:",n_of_cells_x, n_of_cells_y, n_of_cells_z)
    
    cell_size_x, cell_size_y, cell_size_z = find_sizes(n_of_cells_x, n_of_cells_y, n_of_cells_z, geometry_data_object)
    
    precision = 0.1
    #print("CELL SIZE AFTER ADJUSTEMENTS:",(cell_size_x), (cell_size_y), (cell_size_z))
    # if __debug__:
        
    #     for size,quantum in zip([cell_size_x,cell_size_y,cell_size_z],[quantum_x,quantum_y,quantum_z]):
    #         print(abs(size*(1/precision) - quantum),precision)
    #         assert abs(size*(1/precision) - quantum)<=precision
            

    
    n_materials = length(dictData["STLList"])
    
    mesher_output = fill(false, (n_materials,n_of_cells_x, n_of_cells_y, n_of_cells_z))

    mapping_ids_to_materials = Dict()

    counter_stl_files = 1
    for mesh_id in meshes
        
        #@assert meshes[mesh_id] isa Mesh
        #print("voxeling",mesh_id)

        println(mesh_id.first)



        mesher_output[counter_stl_files,:,:,:] = voxelize(n_of_cells_x, n_of_cells_y, n_of_cells_z, meshes[mesh_id.first], geometry_data_object)

        mapping_ids_to_materials[counter_stl_files]=mesh_id
        counter_stl_files+=1
    end

    id_mats_keep=zeros(Int64, n_materials)
    
    id_mats_keep[1]=0
    
    solve_overlapping(n_of_cells_x, n_of_cells_y, n_of_cells_z, n_materials, id_mats_keep, mesher_output)
        



    origin_x = geometry_data_object["meshXmin"]*1e-3
    origin_y = geometry_data_object["meshYmin"]*1e-3
    origin_z = geometry_data_object["meshZmin"]*1e-3


    # assert(isinstance(mesher_output, np.ndarray))
    # @assert cell_size_x isa Float64
    # @assert cell_size_y isa Float64
    # @assert cell_size_z isa Float64
    # @assert origin_x isa Float64
    # @assert origin_y isa Float64
    # @assert origin_z isa Float64

    # Writing to data.json
    json_file = "outputMesher.json"
    
    return dump_json_data(son_file, counter_stl_files, origin_x, origin_y, origin_z, cell_size_x, cell_size_y, cell_size_z,
            n_of_cells_x, n_of_cells_y, n_of_cells_z, mesher_output, mapping_ids_to_materials) 
end

