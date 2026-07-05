# Setting up
import illustris_julia as il
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE, MMFNEXUS, Statistics, Printf, Plots.Measures

basePath = "/net/virgo01/data/users/folkertsma/MScThesis/data/PracticeData/TNG50-4-Dark/output"

println("Loading parameters and setting up PS-DTFE")

# --- Array and Parameters (Marked as CONST for Julia performance) ---
const zz = [0.0, 2.0, 5.0, 10.0]
const ids = [99, 33, 17, 4]
const m = 0.0186754872146347 * 10^10 # DM particle mass in Msun/h
const res = 128                      # resolution
const L = 35.0                       # cMpc/h
const Ni = 270
const Range = range(0, L, length=res)
const depth = 5
const sim_box = SimBox(L, Ni)

# Lagrangian coordinates

# Lagrangian coordinates
println("Loading Lagrangian coordinates...")

# NOTE: Make sure `basePath` is defined somewhere above this!
raw_coords_q = il.snapshot.loadSubset(basePath, 0, "dm", "Coordinates")
id_q = il.snapshot.loadSubset(basePath, 0, "dm", "ParticleIDs")
ordering_q = sortperm(id_q)

# Do all the array operations and assign to const ONCE
const coords_q = (raw_coords_q[:, ordering_q]') ./ 1000

println("Lagrangian coordinates loaded")

# --- Important Stuff ---
function ps_dtfe(snap_id)
    index = findfirst(isequal(snap_id), ids)
    redshift = zz[index]

    folder_name = "RESULTS"
    file_name = "data_z=$redshift.jld2"
    full_path = joinpath(folder_name, file_name)
    
    # Check if the file actually exists before trying to load
    if !isfile(full_path)
        println("Warning: $full_path does not exist. Skipping.")
        return
    end

    raw_coords_x = il.snapshot.loadSubset(basePath, snap_id, "dm", "Coordinates")
    id_x = il.snapshot.loadSubset(basePath, snap_id, "dm", "ParticleIDs")
    raw_vels = il.snapshot.loadSubset(basePath, snap_id, "dm", "Velocities")

    ordering_x = sortperm(id_x)

    coords_x = (raw_coords_x[:, ordering_x]') ./ 1000 # cMpc/h
    vels = raw_vels[:, ordering_x]'                   # km sqrt(a) /s

    ps_dtfe_sb = ps_dtfe_subbox(coords_q, coords_x, vels, m, depth, sim_box; N_target=54)
    estimator = ps_dtfe_sb

    @load full_path density_field_3d norm_dens
    println("Data loaded for z=$redshift")

    coords_arr = [[x, y, z] for x in Range, y in Range, z in Range]
    
    velocity_field_3d = velocity_subbox(coords_arr, estimator)
    stream_field_3d = numberOfStreams_subbox(coords_arr, estimator)
    velocity_field_3d_sum = velocitySum_subbox(coords_arr, estimator)

    # --- SAFE SAVING PROTOCOL ---
    temp_path = joinpath(folder_name, "temp_$file_name")
    
    # 1. Save to a temporary file
    jldsave(temp_path; estimator, density_field_3d, norm_dens, velocity_field_3d, stream_field_3d, velocity_field_3d_sum)
    
    # 2. Overwrite the original ONLY after a successful save
    mv(temp_path, full_path, force=true)
    println("Safely overwrote $full_path")
    
    GC.gc()
end

println("Big calculations about to begin..")

for i in ids
    println("Working on snapshot $i ...")
    ps_dtfe(i)
end