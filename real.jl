## Setting up

import illustris_julia as il
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE, MMFNEXUS, Statistics, Printf, Plots.Measures

basePath = "/net/virgo01/data/users/folkertsma/MScThesis/data/PracticeData/TNG50-4-Dark/output"

println("It's working!")
##
import Pkg
Pkg.resolve()
Pkg.instantiate()qq
## Arrays and Patameters

z = [0.0,0.1,0.2,0.3,0.4,0.5,0.7,1.0,1.5,2.0,3.01,4.01,5.0,6.01,7.01,8.01,9.0,10.0,10.98,11.98,20.05]
ids = [99,91,84,78,72,67,59,50,40,33,25,21,17,13,11,8,6,4,3,2,0]

m = 0.0186754872146347 * 10^10 #DM particle mass in Msun/h
res = 128 #resolution
L = 35.0 #cMpc/h
Ni = 270
Range = range(0, L, length=res)
depth = 5
sim_box = SimBox(L, Ni)

## Deltas arrays

deltas = [] # density contrast for node detection

Omega_0 = 0.3089
Lambda_0 = 0.6911

for zz in z
    Omega = Omega_0 * (1+zz)^3 / (Omega_0 * (1+zz)^3 + Lambda_0)
    d = 18*pi^2 + 82*(Omega - 1) - 39*(Omega - 1)^2
    d = d / Omega
    push!(deltas, d)
end

println(deltas)

#NEXUS Parameters (no delta)
min_node_mass::Real = 1e13 # minimum mass of a node in Msun/h
min_fila_volume::Real = 10 # minimum volume of a filament in (Mpc/h)^3
min_wall_volume::Real = 10 # minimum volume of a wall in (Mpc/h)^3
R0::Real = 0.5 #minimum smoothing scale in Mpc/h
filter_parse = 6 #max n in min_scale*(√2)^n, starting at n=0
level::Symbol = :info # verbose level

## Lagrangian coordinates

coords_q = il.snapshot.loadSubset(basePath, 0, "dm", "Coordinates")
id_q = il.snapshot.loadSubset(basePath, 0, "dm", "ParticleIDs")
ordering_q = sortperm(id_q)
coords_q = coords_q[:,ordering_q]
coords_q = coords_q'
coords_q = coords_q ./ 1000 #convert to cMpc/h
## Important Stuff 
function ps_dtfe(snap_id)

    index = findfirst(isequal(snap_id), ids)
    redshift = z[index]

    folder_name = "RESULTS"
    file_name = "data_z=$redshift"
    full_path = joinpath(folder_name, file_name)

    coords_x = il.snapshot.loadSubset(basePath, snap_id, "dm", "Coordinates")
    id_x = il.snapshot.loadSubset(basePath, snap_id, "dm", "ParticleIDs")
    vels = il.snapshot.loadSubset(basePath, snap_id, "dm", "Velocities")
    ordering_x = sortperm(id_x)
    coords_x = coords_x[:,ordering_x]
    vels = vels[:,ordering_x]
    coords_x = coords_x'
    coords_x = coords_x ./ 1000 #convert to cMpc/h
    vels = vels' #km sqrt(a) /s

    ps_dtfe_sb = ps_dtfe_subbox(coords_q, coords_x, vels, m, depth, sim_box; N_target=54)
    estimator = ps_dtfe_sb
    jldsave(full_path; estimator)

    coords_arr = [[x, y, z] for x in Range, y in Range, z in Range]
    density_field_3d = density_subbox(coords_arr, estimator)
    mean = sum(density_field_3d) / length(density_field_3d)
    norm_dens = density_field_3d / mean #normalized denisty as in 1 + delta
    jldsave(full_path; density_field_3d, norm_dens)
end

function nexus(snap_id)

    index = findfirst(isequal(snap_id), ids)
    redshift = z[index]

    folder_name = "RESULTS"
    file_name = "data_z=$redshift.jld2"
    full_path = joinpath(folder_name, file_name)
    @load full_path norm_dens

    Δ::Real = deltas[index]
    N = res
    totalMass = m * 270^3 # total mass contained in simulation box in Msun/h
    MMF_node, MMF_filament, MMF_wall, MMF_void = NEXUS_Plus(norm_dens, N, L, totalMass;  filter_parse = filter_parse, Δ = Δ, min_node_mass = 1e13, min_fila_volume = min_fila_volume, min_wall_volume = min_wall_volume, R0 = R0, level = level);

    jldsave(full_path; MMF_node, MMF_filament, MMF_wall, MMF_void)
end

## Running the function
ps_dtfe(33)
## Running the other function
nexus(0)
## Saving in HDF5

function save_nexus_results(z_array, filepath="nexus_outputs.h5")
    println("Creating HDF5 file with redshift groups...")
    
    h5open(filepath, "w") do file
        # 1. Global Metadata for the whole simulation
        attributes(file)["algorithms"] = "PS-DTFE and MMF-Nexus"
        attributes(file)["box_size_pixels"] = 128
        attributes(file)["box_size_cMpc/h"] = 35
        
        # 2. Loop through every redshift in your array
        for redshift in z_array
            # Define the path to the JLD2 file saved by your DTFE function
            folder_name = "RESULTS"
            file_name = "data_z=$redshift.jld2" 
            jld2_path = joinpath(folder_name, file_name)
            
            # Safety check: skip this redshift if the file hasn't been generated yet
            if !isfile(jld2_path)
                println("Warning: Could not find $jld2_path. Skipping z=$redshift...")
                continue 
            end
            
            println("Loading and converting data for z = $redshift...")
            
            # Load the variables back into memory from the JLD2 file
            # Note: The variable names must match exactly what you used in jldsave
            if in(redshift, [0.0, 2.0, 5.0, 10.0])
                @load jld2_path density_field_3d norm_dens MMF_node MMF_filament MMF_wall MMF_void velocity_field_3d stream_field_3d velocity_field_3d_sum
            else
                @load jld2_path density_field_3d norm_dens MMF_node MMF_filament MMF_wall MMF_void 
            end
            # 3. Dynamically create the group name 
            group_name = "z_$redshift"
            g = create_group(file, group_name)
            
            # Attach metadata to this specific group
            attributes(g)["redshift"] = Float64(redshift)
            
            # 4. Write the arrays into this newly created group
            write(g, "nodes", Int8.(MMF_node))
            write(g, "filaments", Int8.(MMF_filament))
            write(g, "walls", Int8.(MMF_wall))
            write(g, "voids", Int8.(MMF_void))
            write(g, "normalized_density", norm_dens)
            write(g, "density", density_field_3d)

            if in(redshift, [0.0, 2.0, 5.0, 10.0])
                write(g, "vel_fld", Int8.(velocity_field_3d))
                write(g, "vel_fld_sum", Int8.(velocity_field_3d_sum))
                write(g, "stream_fld", Int8.(stream_field_3d))
            end
        end
    end
    
    println("Successfully saved all automated Nexus data to $filepath")
end

# Run the function, passing in your global 'z' array
save_nexus_results([0.0,0.1,0.2,0.3,0.4,0.5,0.7,1.0,1.5,2.0,3.01,4.01,5.0,6.01,7.01,8.01,9.0,10.0])

## Fixing Nexus

folder_name = "RESULTS"
file_name = "data_z=0.0"
full_path = joinpath(folder_name, file_name)

@load full_path estimator density_field_3d norm_dens MMF_node MMF_filament MMF_wall MMF_void

Δ::Real = deltas[1]
N = res
totalMass = m * 270^3 # total mass contained in simulation box in Msun/h
MMF_node, MMF_filament, MMF_wall, MMF_void = NEXUS_Plus(norm_dens, N, L, totalMass;  filter_parse = filter_parse, Δ = Δ, min_node_mass = min_node_mass, min_fila_volume = min_fila_volume, min_wall_volume = min_wall_volume, R0 = R0, level = level);

folder_name = "RESULTS"
file_name = "data_z=0.0"
full_path = joinpath(folder_name, file_name)
jldsave(full_path; estimator, density_field_3d, norm_dens, MMF_node, MMF_filament, MMF_wall, MMF_void)

##

folder_name = "RESULTS"
file_name = "data_z=0.0.jld2"
full_path = joinpath(folder_name, file_name)

@load full_path estimator density_field_3d norm_dens MMF_node MMF_filament MMF_wall MMF_void

density_field_3d = density_field_3d .* 10^10

jldsave(full_path; estimator, density_field_3d, norm_dens, MMF_node, MMF_filament, MMF_wall, MMF_void)

