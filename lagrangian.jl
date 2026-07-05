## Setting up
import illustris_julia as il
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE, MMFNEXUS, Statistics, Printf, Plots.Measures

println("It's working!")
## Comparison Plot

snap_id = 33 # z=2

basePath = "/net/virgo01/data/users/folkertsma/MScThesis/data/PracticeData/TNG50-4-Dark/output"
header_now = il.groupcat.loadHeader(basePath, 99)

#header_then = il.groupcat.loadHeader(basePath, 33)
#z = header_then["Redshift"]

dm_pos_0 = il.snapshot.loadSubset(basePath,99,"dm",["Coordinates"]);
dm_pos_Mpc_0 = dm_pos_0 ./ 1000 #convert to cMpc/h
dm_pos_2 = il.snapshot.loadSubset(basePath,33,"dm",["Coordinates"]);
dm_pos_Mpc_2 = dm_pos_2 ./ 1000 #convert to cMpc/h

shared_clims = (1, 1000000)

h2 = histogram2d(dm_pos_Mpc_0[2,:], dm_pos_Mpc_0[1,:], colorbar_scale = :log10,
     nbins=500, title="Dark Matter at z = 0",
     xlim=[0,35], ylim=[0,35], xlabel="x [cMpc/h]", ylabel="y [cMpc/h]", margin = 10mm,
     aspect_ratio=:equal, colorbar_title = "Particle Count", clims = shared_clims)


h1 = histogram2d(dm_pos_Mpc_2[2,:], dm_pos_Mpc_2[1,:], colorbar_scale = :log10,
     nbins=500, title="Dark Matter at z = 2", 
     xlim=[0,35], ylim=[0,35], xlabel="x [cMpc/h]", ylabel="y [cMpc/h]", margin = 10mm,
     aspect_ratio=:equal, colorbar_title = "Particle Count", clims = shared_clims)

plot(h1, h2, layout = (1, 2), size = (1200, 500))

## PS_DTFE for z=2

snap_id = 33 # z=2
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
jldsave("estimator_z=2.jld2"; ps_dtfe_sb)

## Recovering the (3D) denisity Field
 
L = 35000.0 #ckpc/h

Range = 0:(L/128):L # 128 pixels across
ps_dtfe_sb = load("estimator.jld2")["ps_dtfe_sb"] #make sure to load correct estimator!

coords_arr = [[x, y, z] for x in Range, y in Range, z in Range]
density_field_3d = density_subbox(coords_arr, ps_dtfe_sb)
jldsave("my_3d_density_field.jld2"; density_field_3d)

#coords_arr = [[x, y, z] for x = L/2, y in Range, z in Range]
#density_field_slice = density_subbox(coords_arr, ps_dtfe_sb)
#jldsave("slice_density_field_128.jld2"; density_field_slice)

## Plotting one slice

L = 35.0 #cMpc/h

Range = 0:(L/128):L # 128 pixels across

slice = load("slice_density_field_128.jld2")["density_field_slice"]
mean = sum(slice) / length(slice)
norm_dens = slice/mean

heatmap(Range, Range, log10.(norm_dens),aspect_ratio=:equal, 
        c=:grays, title="Density Field for z = 0, res = 128", 
        xlabel="[cMpc/h]", ylabel="[cMpc/h]", margin = 10mm,
        xlims = (0,L), ylims = (0,L), colorbar_title = "log₁₀(1+δ)")

## Loading Density Fields

d_f_3d = load("my_3d_density_field.jld2")["density_field_3d"]
mean = sum(d_f_3d) / length(d_f_3d)
norm_dens = d_f_3d/mean #"normalized denisty as in 1 + delta"

d_f_3d_z = load("my_3d_density_field_z=2.jld2")["density_field_3d"]
mean_z = sum(d_f_3d_z) / length(d_f_3d_z)
norm_dens_z = d_f_3d_z/mean_z #"normalized denisty as in 1 + delta"

L = 35.0 #cMpc/h
Range = 0:(L/100):L #cMpc/h

depth = size(norm_dens, 3)
proj = sum(norm_dens, dims=3)[:, :, 1] ./ depth

depth_z = size(norm_dens_z, 3)
proj_z = sum(norm_dens_z, dims=3)[:, :, 1] ./ depth_z

h2 = heatmap(Range, Range, log10.(proj),aspect_ratio=:equal, 
        c=:grays, title="Projected Overdensity Field for z = 0", 
        xlabel="[cMpc/h]", ylabel="[cMpc/h]", margin = 10mm,
        xlims = (0,L), ylims = (0,L), colorbar_title = "log₁₀(1+δ)")

h1 = heatmap(Range, Range, log10.(proj_z),aspect_ratio=:equal, 
        c=:grays, title="Projected Overdensity Field for z = 2", 
        xlabel="[cMpc/h]", ylabel="[cMpc/h]", margin = 10mm,
        xlims = (0,L), ylims = (0,L), colorbar_title = "log₁₀(1+δ)")

plot(h1, h2, layout = (1, 2), size = (1200, 500))

## Nexus implementation
Δ::Real = 370. # density contrast for node detection
min_node_mass::Real = 1e13 # minimum mass of a node in Msun/h
min_fila_volume::Real = 10 # minimum volume of a filament in (Mpc/h)^3
min_wall_volume::Real = 10 # minimum volume of a wall in (Mpc/h)^3
R0::Real = 0.5 #minimum smoothing scale in Mpc/h
filter_parse = 6 #max n in min_scale*(√2)^n, starting at n=0
level::Symbol = :info # verbose level

d_f_3d = load("my_3d_density_field.jld2")["density_field_3d"]
mean = sum(d_f_3d) / length(d_f_3d)
norm_dens = d_f_3d/mean #"normalized denisty as in 1 + delta"

d_f_3d_z = load("my_3d_density_field_z=2.jld2")["density_field_3d"]
mean_z = sum(d_f_3d_z) / length(d_f_3d_z)
norm_dens_z = d_f_3d_z/mean_z #"normalized denisty as in 1 + delta"

N = 101 # number of gridpoints per dimension
L = 35. # Box size in cMpc/h
totalMass = 2.8 * 10^8 * 270^3 # total mass contained in simulation box in Msun/h 

MMF_node_z, MMF_filament_z, MMF_wall_z, MMF_void_z = NEXUS_Plus(norm_dens_z, N, L, totalMass; filter_parse = filter_parse, Δ = Δ, min_node_mass = min_node_mass, min_fila_volume = min_fila_volume, min_wall_volume = min_wall_volume, R0 = R0, level = level);

MMF_node, MMF_filament, MMF_wall, MMF_void = NEXUS_Plus(norm_dens, N, L, totalMass;  filter_parse = filter_parse, Δ = Δ, min_node_mass = min_node_mass, min_fila_volume = min_fila_volume, min_wall_volume = min_wall_volume, R0 = R0, level = level);


## Animation Nexus
N = size(norm_dens_z, 1)
L = 35.0           
spatial_coords = range(0, stop=L, length=N)
mini = minimum(norm_dens_z)
maxi= maximum(norm_dens_z)
anim = @animate for slice_index in 1:size(norm_dens_z,3)
        heatmap(spatial_coords,spatial_coords,log10.(norm_dens_z[:, :, slice_index]), 
                aspect_ratio=:equal, c=:grays, 
                title=@sprintf("Density Field at z=2 for slice %d", slice_index), 
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

mp4(anim, "nexus_cinema_z=2.mp4", fps = 15)

## Plotting individual features

voids = norm_dens .* MMF_void
filaments = norm_dens .* MMF_filament
nodes = norm_dens .* MMF_node
walls = norm_dens .* MMF_wall

voids_z = norm_dens_z .* MMF_void_z
filaments_z = norm_dens_z .* MMF_filament_z
nodes_z = norm_dens_z .* MMF_node_z
walls_z = norm_dens_z .* MMF_wall_z

slice_id = 50

h1 = heatmap(log10.(voids[:,:,slice_id]), title = "Voids")
h2 = heatmap(log10.(walls[:,:,slice_id]), title = "Walls")
h3 = heatmap(log10.(nodes[:,:,slice_id]), title = "Nodes")
h4 = heatmap(log10.(filaments[:,:,slice_id]), title = "Filaments")

h5 = heatmap(log10.(voids_z[:,:,slice_id]), title = "Voids")
h6 = heatmap(log10.(walls_z[:,:,slice_id]), title = "Walls")
h7 = heatmap(log10.(nodes_z[:,:,slice_id]), title = "Nodes")
h8 = heatmap(log10.(filaments_z[:,:,slice_id]), title = "Filaments")

#plot(h1, h2, h3, h4, layout = (2, 2), plot_title = "z=0")

plot(h5, h6, h7, h8, layout = (2, 2), plot_title = "z=2")

## Statistics z=0

voids = norm_dens .* MMF_void
filaments = norm_dens .* MMF_filament
nodes = norm_dens .* MMF_node
walls = norm_dens .* MMF_wall

voids_z = norm_dens_z .* MMF_void_z
filaments_z = norm_dens_z .* MMF_filament_z
nodes_z = norm_dens_z .* MMF_node_z
walls_z = norm_dens_z .* MMF_wall_z

histogram(log10.(vec(voids)), label = "Voids", title = "z=0", normalize = :pdf, xlabel="log₁₀(1+δ)", ylabel="PDF")
histogram!(log10.(vec(filaments)), label = "Filaments", normalize = :pdf)
histogram!(log10.(vec(walls)), label = "Walls", normalize = :pdf)
histogram!(log10.(vec(nodes)), label = "Nodes", normalize = :pdf)

## Statistics z=2

histogram(log10.(vec(voids_z)), label = "Voids", title = "z=2",normalize = :pdf, xlabel="log₁₀(1+δ)", ylabel="PDF")
histogram!(log10.(vec(filaments_z)), label = "Filaments", normalize = :pdf)
histogram!(log10.(vec(walls_z)), label = "Walls", normalize = :pdf)
histogram!(log10.(vec(nodes_z)), label = "Nodes", normalize = :pdf)

## Files for Python

function save_nexus_results(filepath="nexus_outputs.h5")
    println("Creating HDF5 file with redshift groups...")
    
    h5open(filepath, "w") do file
        # 1. Global Metadata for the whole simulation
        attributes(file)["algorithm"] = "MMF/Nexus"
        attributes(file)["box_size_mpc"] = 101.0 # Adjust to your sim
        
        # ==========================================
        # 2. Create a "Folder" (Group) for z = 0
        # ==========================================
        g_z0 = create_group(file, "z_0")
        attributes(g_z0)["redshift"] = 0.0 # Attach metadata just to this folder
        
        # Save the 4 arrays inside the z=0 group
        write(g_z0, "nodes", nodes)
        write(g_z0, "filaments", filaments)
        write(g_z0, "walls", walls)
        write(g_z0, "voids", voids)
        write(g_z0, "density", norm_dens)
        
        # ==========================================
        # 3. Create a "Folder" (Group) for z = 2
        # ==========================================
        g_z2 = create_group(file, "z_2")
        attributes(g_z2)["redshift"] = 2.0
        
        # Save the 4 arrays inside the z=2 group
        write(g_z2, "nodes", nodes_z)
        write(g_z2, "filaments", filaments_z)
        write(g_z2, "walls", walls_z)
        write(g_z2, "voids", voids_z)
        write(g_z2, "density", norm_dens_z)
    end
    
    println("Successfully saved Nexus data to $filepath")
end

# Run the function to save
save_nexus_results()

## Animation

N = size(norm_dens, 1)
L = 35.0           
spatial_coords = range(0, stop=L, length=N)
mini = minimum(norm_dens_z)
maxi= maximum(norm_dens_z)

anim = @animate for slice_index in 1:size(norm_dens,3)
        h1 = heatmap(log10.(voids_z[:,:,slice_index]), title = "Voids")
        h2 = heatmap(log10.(walls_z[:,:,slice_index]), title = "Walls")
        h3 = heatmap(log10.(nodes_z[:,:,slice_index]), title = "Nodes")
        h4 = heatmap(log10.(filaments_z[:,:,slice_index]), title = "Filaments")
        plot(h1, h2, h3, h4, layout = (2, 2), plot_title = @sprintf("z=2 for slice %d", slice_index))

end

mp4(anim, "stuff_at_z=2.mp4", fps = 15)

# CDF
# 240 for delta! (for specific z)
# R0 at 0.5
# filter_parse = 6

## Bla

k = load("my_3d_density_field.jld2")["density_field_3d"]
size(k)

## Testing out different Nexus
Δ::Real = 331. # density contrast for node detection
min_node_mass::Real = 1e13 # minimum mass of a node in Msun/h
min_fila_volume::Real = 10 # minimum volume of a filament in (Mpc/h)^3
min_wall_volume::Real = 10 # minimum volume of a wall in (Mpc/h)^3
R0::Real = 0.5 #minimum smoothing scale in Mpc/h
filter_parse = 6 #max n in min_scale*(√2)^n, starting at n=0
level::Symbol = :info # verbose level

folder_name = "RESULTS"
file_name = "data_z=0.0"
full_path = joinpath(folder_name, file_name)

norm_dens = load(full_path)["norm_dens"]
m = 0.0186754872146347 * 10^10 #DM particle mass in Msun/h

N = 128 # number of gridpoints per dimension
L = 35.0 # Box size in cMpc/h
totalMass = m * 270^3 # total mass contained in simulation box in Msun/h

MMF_node, MMF_filament, MMF_wall, MMF_void = NEXUS_Plus(norm_dens, N, L, totalMass;  filter_parse = filter_parse, Δ = Δ, min_node_mass = min_node_mass, min_fila_volume = min_fila_volume, min_wall_volume = min_wall_volume, R0 = R0, level = level);

N = size(norm_dens, 1)
L = 35.0           
spatial_coords = range(0, stop=L, length=N)
mini = minimum(norm_dens)
maxi= maximum(norm_dens)
anim = @animate for slice_index in 1:size(norm_dens,3)
        heatmap(spatial_coords,spatial_coords,log10.(norm_dens[:, :, slice_index]), 
                aspect_ratio=:equal, c=:grays, 
                title=@sprintf("Density Field at z=0 for slice %d", slice_index), 
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

mp4(anim, "nexus_test.mp4", fps = 15)

##
voids = norm_dens .* MMF_void
filaments = norm_dens .* MMF_filament
nodes = norm_dens .* MMF_node
walls = norm_dens .* MMF_wall

histogram(log10.(vec(voids)), label = "Voids", title = "z=0", normalize = :pdf, xlabel="log₁₀(1+δ)", ylabel="PDF")
histogram!(log10.(vec(filaments)), label = "Filaments", normalize = :pdf)
histogram!(log10.(vec(walls)), label = "Walls", normalize = :pdf)
histogram!(log10.(vec(nodes)), label = "Nodes", normalize = :pdf)

## opening density
redshift = 1.5

folder_name = "RESULTS"
file_name = "data_z=$redshift.jld2" 
jld2_path = joinpath(folder_name, file_name)
@load jld2_path density_field_3d norm_dens 

N = size(norm_dens, 1)
L = 35.0           
spatial_coords = range(0, stop=L, length=N)
mini = minimum(norm_dens)
maxi= maximum(norm_dens)

slice_index = 64

anim = @animate for slice_index in 1:size(norm_dens,3)
        heatmap(spatial_coords,spatial_coords,log10.(norm_dens[:, :, slice_index]), 
        aspect_ratio=:equal, c=:grays, 
        title="Density Field at z=$redshift for slice $slice_index", 
        xlabel="x [Mpc/h]", ylabel="y [Mpc/h]",xlims=(0,L), ylims=(0,L),
        colorbar_title="log₁₀(1+δ)", 
        clims = (log10(mini), log10(maxi)))
end

mp4(anim, "density_z=$redshift.mp4", fps = 15)
## 
redshift = 2.0
estimator = load("estimator_z=2.jld2")["ps_dtfe_sb"]

folder_name = "RESULTS"
file_name = "data_z=$redshift.jld2"
full_path = joinpath(folder_name, file_name)

res = 128
L = 35.0 #cMpc/h
Ni = 270
Range = range(0, L, length=res)

coords_arr = [[x, y, z] for x in Range, y in Range, z in Range]
density_field_3d = density_subbox(coords_arr, estimator)
mean = sum(density_field_3d) / length(density_field_3d)
norm_dens = density_field_3d / mean #normalized denisty as in 1 + delta

jldsave(full_path; estimator, density_field_3d, norm_dens)

## density evolution

zz = [0.0,0.1,0.2,0.3,0.4,0.5,0.7,1.0,1.5,2.0,3.01,4.01,5.0,6.01,7.01,8.01,9.0,10.0,10.98,11.98,20.05]

folder_name = "RESULTS"

N = 128
L = 35.0           
spatial_coords = range(0, stop=L, length=N)

slice_index = 64

slices = []
for z in zz
        file_name = "data_z=$z.jld2" 
        jld2_path = joinpath(folder_name, file_name)
        @load jld2_path norm_dens
        push!(slices, norm_dens[:,:,slice_index])
end

slices = stack(slices)

mini = minimum(slices)
maxi = maximum(slices)

anim = @animate for i in range(1,21)
        z = size(slices)
        file_name = "data_z=$z.jld2" 
        jld2_path = joinpath(folder_name, file_name)
        @load jld2_path norm_dens
        heatmap(spatial_coords,spatial_coords,log10.(slices[i]), 
        aspect_ratio=:equal, c=:grays, 
        title="Density Field at z=$z for slice $slice_index", 
        xlabel="x [Mpc/h]", ylabel="y [Mpc/h]",xlims=(0,L), ylims=(0,L),
        colorbar_title="log₁₀(1+δ)", 
        clims = (log10(mini), log10(maxi))
        )
end

mp4(anim, "density_evolution.mp4", fps = 5)

##
using JLD2, Plots # Making sure your dependencies are active

zz = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.01, 4.01, 5.0, 6.01, 7.01, 8.01, 9.0, 10.0, 10.98, 11.98, 20.05]
zz_reversed = reverse(zz)

folder_name = "RESULTS"

N = 128
L = 35.0           
spatial_coords = range(0, stop=L, length=N)

slice_index = 64

slices_vec = []
for z in zz
        file_name = "data_z=$z.jld2" 
        jld2_path = joinpath(folder_name, file_name)
        @load jld2_path norm_dens
        push!(slices_vec, norm_dens[:, :, slice_index])
end
slices = stack(slices_vec)
slices_reversed = reverse(slices)

mini = minimum(slices)
maxi = maximum(slices)

anim = @animate for i in 1:length(zz)
        z = zz_reversed[i] 
        
        current_slice = slices_reversed[:, :, i]
        
        heatmap(spatial_coords, spatial_coords, log10.(current_slice), 
            aspect_ratio=:equal, 
            c=:grays, 
            title="Density Field at z=$z for slice $slice_index", 
            xlabel="x [Mpc/h]", 
            ylabel="y [Mpc/h]",
            xlims=(0, L), 
            ylims=(0, L),
            colorbar_title="log₁₀(1+δ)", 
            clims=(log10(mini), log10(maxi))
        )
end

mp4(anim, "density_evolution.mp4", fps = 5)