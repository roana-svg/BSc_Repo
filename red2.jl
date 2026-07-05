import illustris_julia as il
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE, MMFNEXUS, Statistics, Printf, Plots.Measures

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


## ufurf
import illustris_julia as il
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE, MMFNEXUS, Statistics, Printf, Plots.Measures


data = load("RESULTS/data_z=20.05.jld2")

for key in keys(data)
    # Get the size in bytes
    size_bytes = Base.summarysize(data[key])
    
    # Convert to MB for easier reading
    size_mb = size_bytes / (1024 * 1024)
    
    # Print the result rounded to 2 decimal places
    println("Key: $key | Size: $(round(size_mb, digits=2)) MB")
end

##

@load "RESULTS/data_z=20.05.jld2" estimator

# Look at all the accessible field names inside the estimator
println("Fields inside the estimator:")
println(propertynames(estimator))

# Print the entire internal structure and values (this will output a lot of text!)
dump(estimator)

##
strr = load("/Users/users/roana/roana/BSc_Thesis/RESULTS/data_z=0.0.jld2")["stream_field_3d"]

L = 35.0 #cMpc/h
Range = 0:(L/128):L #cMpc/h


#depth_z = size(norm_dens_z, 3)
#proj_z = sum(norm_dens_z, dims=3)[:, :, 1] ./ depth_z

heatmap(Range, Range, log10.(strr[:,:,67]),aspect_ratio=:equal, 
        c=:inferno, title="Stream Field for z = 20.05, slice 67", 
        xlabel="[cMpc/h]", ylabel="[cMpc/h]", margin = 10mm,
        xlims = (0,L), ylims = (0,L), colorbar_title = "log(#streams)")