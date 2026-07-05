## Setting up
import illustris_julia as il
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE, MMFNEXUS, Statistics, Printf, Plots.Measures

println("It's working!")

##
stream = load("/Users/users/roana/roana/BSc_Thesis/RESULTS/data_z=0.0.jld2")["stream_field_3d"]

L = 35.0 #cMpc/h
Range = 0:(L/128):L #cMpc/h


#depth_z = size(norm_dens_z, 3)
#proj_z = sum(norm_dens_z, dims=3)[:, :, 1] ./ depth_z

heatmap(Range, Range, log10.(stream[:,:,67]),aspect_ratio=:equal, 
        c=:inferno, title="Stream Field for z = 20.05, slice 67", 
        xlabel="[cMpc/h]", ylabel="[cMpc/h]", margin = 10mm,
        xlims = (0,L), ylims = (0,L), colorbar_title = "log(#streams)")

##

vels = load("/Users/users/roana/roana/BSc_Thesis/RESULTS/data_z=20.05.jld2")["velocity_field_3d"]
n_d = load("/Users/users/roana/roana/BSc_Thesis/RESULTS/data_z=20.05.jld2")["norm_dens"]

Vx = vels[:,:,:,1]
Vy = vels[:,:,:,2]

L = 35.0 # cMpc/h
Nx, Ny, Nz = size(Vx)
k = 67

# 1. Create a spatial range of exactly 128 points
Range = range(0, L, length=Nx)

# 2. Plot the base heatmap
p = heatmap(Range, Range, log10.(n_d[:,:,k]), 
        aspect_ratio=:equal, 
        c=:inferno, 
        title="Stream & Velocity Field (z = 0, slice $k)", 
        xlabel="[cMpc/h]", ylabel="[cMpc/h]", margin=10mm,
        xlims=(0,L), ylims=(0,L), colorbar_title="log(1+δ)")

# 3. Define a stride and get physical coordinates for the quiver
stride = 4 
x_range = 1:stride:Nx
y_range = 1:stride:Ny

# Map indices to physical coordinates for the quiver plot
x_coords = [Range[x] for x in x_range, y in y_range] |> vec
y_coords = [Range[y] for x in x_range, y in y_range] |> vec

# 4. Slice, flatten, AND SCALE the velocity components
scale_factor = 0.02 # <-- TWEAK THIS: smaller number = shorter arrows
u_slice = (Vx[x_range, y_range, k] .* scale_factor) |> vec
v_slice = (Vy[x_range, y_range, k] .* scale_factor) |> vec

# 4. Calculate magnitudes safely
magnitudes = sqrt.(u_slice.^2 .+ v_slice.^2)

# DEFENSE 1: Purge any NaNs from the data (replace them with 0.0 so they don't break the math)
replace!(magnitudes, NaN => 0.0)
replace!(u_slice, NaN => 0.0)
replace!(v_slice, NaN => 0.0)

# 5. Normalize vectors
u_norm = u_slice ./ (magnitudes .+ 1e-6) # 1e-6 is safer than eps() for Float32 simulation data
v_norm = v_slice ./ (magnitudes .+ 1e-6)

# Scale them
fixed_length = 0.8
u_plot = u_norm .* fixed_length
v_plot = v_norm .* fixed_length

# 6. Map magnitudes to Opacity (Alpha) Safely
max_mag = maximum(magnitudes)

# If max_mag is exactly 0 (empty slice), avoid division by zero
if max_mag == 0.0
    max_mag = 1.0 
end

alphas = 0.1 .+ 0.9 .* (magnitudes ./ max_mag)
arrow_colors = [RGBA(1.0, 1.0, 1.0, a) for a in alphas]

# DEFENSE 2: The row-vector trick for linecolor using `reshape`
quiver!(p, x_coords, y_coords, quiver=(u_plot, v_plot), 
        linecolor=reshape(arrow_colors, 1, :), lw=1.5)

display(p)

##

cool = load("/Users/users/spirov/Blk/JuliaDTFE/saves/2D3200.jld2")

heatmap(log10.(cool), 
        aspect_ratio=:equal, 
        c=:inferno, 
        title="HEHE", 
        xlabel="[cMpc/h]", ylabel="[cMpc/h]", margin=10mm,
        xlims=(0,L), ylims=(0,L), colorbar_title="log(1+δ)")