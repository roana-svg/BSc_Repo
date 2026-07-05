# Setting up

import illustris_julia as il
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE, MMFNEXUS, Statistics, Printf, Plots.Measures

basePath = "/net/virgo01/data/users/folkertsma/MScThesis/data/PracticeData/TNG50-4-Dark/output"

println("Loading parameters and setting up PS-DTFE")

# Arrays and Patameters

#zz = [0.0,0.1,0.2,0.3,0.4,0.5,0.7,1.0,1.5,2.0,3.01,4.01,5.0,6.01,7.01,8.01,9.0,10.0,10.98,11.98,20.05]
#ids = [99,91,84,78,72,67,59,50,40,33,25,21,17,13,11,8,6,4,3,2,0]

#zz = [0.1,0.2,0.3,0.4,0.5,0.7,1.0,1.5,2.0,3.01,4.01,5.0,6.01,7.01,8.01,9.0,10.0,10.98,11.98]
#ids = [91,84,78,72,67,59,50,40,33,25,21,17,13,11,8,6,4,3,2]

ids = [0]
zz = [20.05]

m = 0.0186754872146347 * 10^10 #DM particle mass in Msun/h
res = 128 #resolution
L = 35.0 #cMpc/h
Ni = 270
Range = range(0, L, length=res)
depth = 5
sim_box = SimBox(L, Ni)

# Lagrangian coordinates

coords_q = il.snapshot.loadSubset(basePath, 0, "dm", "Coordinates")
id_q = il.snapshot.loadSubset(basePath, 0, "dm", "ParticleIDs")
ordering_q = sortperm(id_q)
coords_q = coords_q[:,ordering_q]
coords_q = coords_q'
coords_q = coords_q ./ 1000

println("Lagrangian coordinates loaded")

# Important Stuff 
function ps_dtfe(snap_id)

    index = findfirst(isequal(snap_id), ids)
    redshift = zz[index]

    folder_name = "RESULTS"
    file_name = "data_z=$redshift.jld2"
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

    coords_arr = [[x, y, z] for x in Range, y in Range, z in Range]
    density_field_3d = density_subbox(coords_arr, estimator)
    meann = sum(density_field_3d) / length(density_field_3d)
    norm_dens = density_field_3d / meann #normalized denisty as in 1 + delta

    jldsave(full_path; estimator, density_field_3d, norm_dens)
    
    GC.gc()
end

println("Big calculations about to begin")

for i in ids
    println("Working on snapshot $i ...")
    ps_dtfe(i)
end