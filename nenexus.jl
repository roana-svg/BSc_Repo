import illustris_julia as il
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE, MMFNEXUS, Statistics, Printf, Plots.Measures

println("Loading parameters and setting up NEXUS")

#const z = [0.0,0.1,0.2,0.3,0.4,0.5,0.7,1.0,1.5,2.0,3.01,4.01,5.0,6.01,7.01,8.01,9.0,10.0,10.98,11.98,20.05]
#const ids = [99,91,84,78,72,67,59,50,40,33,25,21,17,13,11,8,6,4,3,2,0]

const z = [0.0, 2.0, 5.0, 10.0]
const ids = [99, 33, 17, 4]

const deltas = Float64[] # density contrast for node detection

const Omega_0 = 0.3089
const Lambda_0 = 0.6911

for zz in z
    Omega = Omega_0 * (1+zz)^3 / (Omega_0 * (1+zz)^3 + Lambda_0)
    d = 18*pi^2 + 82*(Omega - 1) - 39*(Omega - 1)^2
    d = d / Omega
    push!(deltas, d)
end

println(deltas)

const m = 0.0186754872146347 * 10^10 #DM particle mass in Msun/h
const res = 128 #resolution
const L = 35.0 #cMpc/h
const Ni = 270
const totalMass = m * Ni^3 # total mass contained in simulation box in Msun/h

const min_node_mass = 1e13 # minimum mass of a node in Msun/h
const min_fila_volume = 10 # minimum volume of a filament in (Mpc/h)^3
const min_wall_volume = 10 # minimum volume of a wall in (Mpc/h)^3
const R0 = 0.5 #minimum smoothing scale in Mpc/h
const filter_parse = 8 #max n in min_scale*(√2)^n, starting at n=0
const level = :info # verbose level

function nexus(snap_id)

    index = findfirst(isequal(snap_id), ids)
    redshift = z[index]

    folder_name = "RESULTS"
    file_name = "data_z=$redshift.jld2"
    full_path = joinpath(folder_name, file_name)
    @load full_path estimator density_field_3d norm_dens

    Δ = deltas[index]

    if redshift > 2
        MMF_node, MMF_filament, MMF_wall, MMF_void = NEXUS_Plus(norm_dens, res, L, totalMass;  filter_parse = filter_parse, Δ = Δ, min_node_mass = 1e12, min_fila_volume = min_fila_volume, min_wall_volume = min_wall_volume, R0 = R0, level = level);
    else
        MMF_node, MMF_filament, MMF_wall, MMF_void = NEXUS_Plus(norm_dens, res, L, totalMass;  filter_parse = filter_parse, Δ = Δ, min_node_mass = min_node_mass, min_fila_volume = min_fila_volume, min_wall_volume = min_wall_volume, R0 = R0, level = level);
    end

    jldsave(full_path;  estimator, density_field_3d, norm_dens, MMF_node, MMF_filament, MMF_wall, MMF_void)
end

const failed_snaps = Int[] 

for i in ids
    println("----------------------------------------")
    println("NEXUSING SNAPSHOT $i ...")
    
    try
        @time nexus(i)
        
    catch e
        println("\n❌ ERROR: Snapshot $i failed!")
        
        showerror(stdout, e, catch_backtrace()) 
        println("\nSkipping to the next snapshot...")
        
        push!(failed_snaps, i)
    end
end

println("========================================")
if isempty(failed_snaps)
    println("✅ SUCCESS: All snapshots processed without errors.")
else
    println("⚠️ WARNING: The run completed, but the following snapshots failed:")
    println(failed_snaps)
    println("Check the log above for specific error traces.")
end
println("========================================")