## Setting up 
import illustris_julia as il
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE

basePath = "/net/virgo01/data/users/folkertsma/MScThesis/data/PracticeData/TNG50-4-Dark/output"

header = il.groupcat.loadHeader(basePath, 99)
#println("Box Size: ", header["BoxSize"])
println(header)

## Plotting
using Plots
dm_pos = il.snapshot.loadSubset(basePath,99,"dm",["Coordinates"]);
histogram2d(dm_pos[1,:], dm_pos[2,:], colorbar_scale = :log10, nbins=500)
histogram2d!(xlim=[0,30000], ylim=[0,30000], xlabel="x [ckpc/h]", ylabel="y [ckpc/h]", aspect_ratio=:equal)


## Good PS-DTFE

using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE

snap_id = 99
m = 0.0186754872146347 # 10^10 Msun/h 
#m = m*(10^10) #Msun/h 

coords_q = il.snapshot.loadSubset(basePath, 0, "dm", "Coordinates")
id_q = il.snapshot.loadSubset(basePath, 0, "dm", "ParticleIDs")
ordering_q = sortperm(id_q)
coords_q = coords_q[:,ordering_q]
coords_q = coords_q'


coords_x = il.snapshot.loadSubset(basePath, snap_id, "dm", "Coordinates")
id_x = il.snapshot.loadSubset(basePath, snap_id, "dm", "ParticleIDs")
vels = il.snapshot.loadSubset(basePath, snap_id, "dm", "Velocities")
ordering_x = sortperm(id_x)
coords_x = coords_x[:,ordering_x]
vels = vels[:,ordering_x]
coords_x = coords_x'
vels = vels'



L = 10000.0  # in ckpc/h
box_min = 0.0
box_max = L

mask = (coords_x[:, 1] .>= box_min) .& (coords_x[:, 1] .<= box_max) .&
       (coords_x[:, 2] .>= box_min) .& (coords_x[:, 2] .<= box_max) .&
       (coords_x[:, 3] .>= box_min) .& (coords_x[:, 3] .<= box_max)

coords_x = coords_x[mask, :]
coords_q = coords_q[mask, :]
vels     = vels[mask, :]

N = size(coords_x, 1)

# Format: [xmin xmax; ymin ymax; zmin zmax], NO SIMBOX
box = [0.0 L; 
       0.0 L; 
       0.0 L]
depth = 5

ps_dtfe = PS_DTFE(coords_q, coords_x, vels, m, depth, box)

Range = 0:(L/200):L # 200 pixels across
density_field = [PhaseSpaceDTFE.density([L/2., y, z], ps_dtfe) for y in Range, z in Range]
## Plotting
heatmap(Range, Range, log10.(density_field), aspect_ratio=:equal, 
        xlims=(0, L), ylims=(0, L), c=:grays, 
        title="Test Patch Density", xlabel="[ckpc/h]", ylabel="[ckpc/h]")

## Subbox PS-DTFE
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE

snap_id = 99
m = 0.0186754872146347 # 10^10 Msun/h 
#m = m*(10^10) #Msun/h 


coords_q = il.snapshot.loadSubset(basePath, 0, "dm", "Coordinates")
id_q = il.snapshot.loadSubset(basePath, 0, "dm", "ParticleIDs")
ordering_q = sortperm(id_q)
coords_q = coords_q[:,ordering_q]
coords_q = coords_q'

coords_x = il.snapshot.loadSubset(basePath, snap_id, "dm", "Coordinates")
id_x = il.snapshot.loadSubset(basePath, snap_id, "dm", "ParticleIDs")
vels = il.snapshot.loadSubset(basePath, snap_id, "dm", "Velocities")
ordering_x = sortperm(id_x)
coords_x = coords_x[:,ordering_x]
vels = vels[:,ordering_x]
coords_x = coords_x'
vels = vels'


L = 35000.0  # in ckpc/h
Ni = 270 
sim_box = SimBox(L, Ni)

depth = 5

ps_dtfe_sb = ps_dtfe_subbox(coords_q, coords_x, vels, m, depth, sim_box; N_target=54)
jldsave("estimator.jld2"; ps_dtfe_sb)
## Density slice
L = 35000.0 
Range = 0:(L/400):L # 400 pixels across
ps_dtfe_sb = load("estimator.jld2")["ps_dtfe_sb"]

coords_arr = [[L/2., y, z] for y in Range, z in Range]
density_field = density_subbox(coords_arr, ps_dtfe_sb)
jldsave("my_density_field.jld2"; density_field)
## Plot Plot
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE
L = 35000.0
Range = 0:(L/400):L
Range_cMpc = Range ./ 1000.0
d_f = load("my_density_field.jld2")["density_field"]
#d_f_3d = load("my_3d_density_field.jld2")["density_field_3d"]

heatmap(Range_cMpc, Range_cMpc, log10.(d_f), 
        aspect_ratio=:equal, 
        c=:grays, 
        title="Density Field", 
        xlabel="[cMpc/h]", 
        ylabel="[cMpc/h]",
        xlims = (0,L/1000),
        ylims = (0,L/1000))
## Density box
Range = 0:(L/100):L # 100 pixels across
#ps_dtfe_sb = load("estimator.jld2")["ps_dtfe_sb"]

coords_arr = [[x, y, z] for x in Range, y in Range, z in Range]
density_field_3d = density_subbox(coords_arr, ps_dtfe_sb)
jldsave("my_3d_density_field.jld2"; density_field_3d)
## Check stuff
d_f_3d = load("my_3d_density_field.jld2")["density_field_3d"]
size(d_f_3d[1,:,:])
##

d_f_3d = load("my_3d_density_field.jld2")["density_field_3d"]

mean = sum(d_f_3d) / length(d_f_3d)
delta = d_f_3d/mean
## Animate
L = 35000.0
Range = 0:(L/100):L
Range_cMpc = Range ./ 1000.0
using Printf
mini = minimum(delta)
maxi= maximum(delta)
anim = @animate for k in 1:size(delta,3)
        heatmap(Range_cMpc, Range_cMpc, log10.(delta[k,:,:]), 
        aspect_ratio=:equal, 
        c=:viridis, 
        title=@sprintf("Density Field for slice %d", k), 
        xlabel="[cMpc/h]", 
        ylabel="[cMpc/h]",
        colorbar_title="log₁₀(1+δ)",
        xlims = (0,L/1000),
        ylims = (0,L/1000),
        clims = (log10(mini), log10(maxi)))
end

mp4(anim, "density_slices_video.mp4", fps = 10)

## Nexus of sorts
using MMFNEXUS, Statistics

N = 101 # number of gridpoints per dimension
L = 35. # Box size in cMpc/h
totalMass = 2.8 * 10^8 * 270^3 # total mass contained in simulation box in Msun/h 

MMF_node, MMF_filament, MMF_wall, MMF_void = NEXUS_Plus(delta, N, L, totalMass);
## Nexus plot
N = size(delta, 1)
L = 35.0           
spatial_coords = range(0, stop=L, length=N)

slice_index = 75

heatmap(spatial_coords,spatial_coords,log10.(delta[:, :, slice_index]), 
        aspect_ratio=:equal, c=:grays, title=@sprintf("Density Field for slice %d", slice_index), 
        xlabel="x [Mpc/h]", ylabel="y [Mpc/h]",xlims=(0,L), 
        colorbar_title="log₁₀(1+δ)", 
        legend = :outerbottom,
        legend_columns = 3,       
        bottom_margin = 15Plots.mm)
contour!(spatial_coords,spatial_coords,MMF_wall[:, :, slice_index], levels=[0.5], color=:green, linewidth=2, label="Walls")
contour!(spatial_coords,spatial_coords,MMF_filament[:, :, slice_index], levels=[0.5], color=:blue, linewidth=2, label="Filaments")
contour!(spatial_coords,spatial_coords,MMF_node[:, :, slice_index], levels=[0.5], color=:red, linewidth=2, label="Nodes")
plot!([NaN], [NaN], color=:green, linewidth=2, label="Walls")
plot!([NaN], [NaN], color=:blue, linewidth=2, label="Filaments")
plot!([NaN], [NaN], color=:red, linewidth=2, label="Nodes")
## Animation Nexus
N = size(delta, 1)
L = 35.0           
spatial_coords = range(0, stop=L, length=N)
using Printf
mini = minimum(delta)
maxi= maximum(delta)
anim = @animate for slice_index in 1:size(delta,3)
        heatmap(spatial_coords,spatial_coords,log10.(delta[:, :, slice_index]), 
                aspect_ratio=:equal, c=:grays, 
                title=@sprintf("Density Field for slice %d", slice_index), 
                xlabel="x [Mpc/h]", ylabel="y [Mpc/h]",xlims=(0,L), ylims=(0,L),
                colorbar_title="log₁₀(1+δ)", 
                clims = (log10(mini), log10(maxi)),
                legend = :outerbottom,
                legend_columns = 3,       
                bottom_margin = 10Plots.mm)
        contour!(spatial_coords,spatial_coords,MMF_wall[:, :, slice_index], levels=[0.5], color=:green, linewidth=2, label="Walls")
        contour!(spatial_coords,spatial_coords,MMF_filament[:, :, slice_index], levels=[0.5], color=:blue, linewidth=2, label="Filaments")
        contour!(spatial_coords,spatial_coords,MMF_node[:, :, slice_index], levels=[0.5], color=:red, linewidth=2, label="Nodes")

        plot!([NaN], [NaN], color=:green, linewidth=2, label="Walls")
        plot!([NaN], [NaN], color=:blue, linewidth=2, label="Filaments   ")
        plot!([NaN], [NaN], color=:red, linewidth=2, label="Nodes")
end

mp4(anim, "nexus_cinema.mp4", fps = 15)
##
using Pkg
Pkg.activate(".")    # This tells Julia to look at the Project.toml in your current folder
Pkg.instantiate()    # This tells Julia to download and install everything listed in it
