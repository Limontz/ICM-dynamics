using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using DataFrames

index = range(1,8,step=1) |> collect
name = ["IT90_0","IT90_0", "IT90_1", "IT90_2", "IT90_3", "IT90_4", "IT92_0", "IT92_1", "IT92_2"]


main = "/home/marco/Scrivania/Ettori_project/"

file1=string(main,"Plots/2D_x_radial_filtering_fluctuation_radial_profile_slope_distribution.txt")
file2=string(main,"Plots/2D_y_radial_filtering_fluctuation_radial_profile_slope_distribution.txt")
file3=string(main,"Plots/2D_z_radial_filtering_fluctuation_radial_profile_slope_distribution.txt")
file4=string(main,"Plots/3D_radial_filtering_fluctuation_radial_profile_slope_distribution.txt")

data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
data3=readdlm(file3,comments=true)
data4=readdlm(file3,comments=true)


density = DataFrame( x=data1[:,1], y = data2[:,1], z=data3[:,1])

temperature = DataFrame( x=data1[:,2], y = data2[:,2], z=data3[:,2])

velocity = DataFrame( x=data1[:,3], y = data2[:,3], z=data3[:,3])

density_3D = data4[:,1]
temperature_3D = data4[:,2]
velocity_3D = data4[:,3]
#println(density)
#println(density[:,1])

# Calculating the average among the 3 direction x,y,z
density_mean = (mean(Array(density[:,:]), dims=2))
density_max, density_min = (maximum(Array(density[:,:]), dims=2)), (minimum(Array(density[:,:]), dims=2))

temperature_mean = (mean(Array(temperature[:,:]), dims=2))
temperature_max, temperature_min = (maximum(Array(temperature[:,:]), dims=2)), (minimum(Array(temperature[:,:]), dims=2))

velocity_mean = (mean(Array(velocity[:,:]), dims=2))
velocity_max, velocity_min = (maximum(Array(velocity[:,:]), dims=2)), (minimum(Array(velocity[:,:]), dims=2))

#indices of perturbed clusters
perturbed_indices = [1,4,5,7,8]
relaxed_indices = [2,3,6]

println("3D density slope (t,p,r) = ", mean(Array(density_3D)), "+-", std(Array(density_3D)),  " ", mean(Array(density_3D[perturbed_indices,:])),  "+-", std(Array(density_3D[perturbed_indices,:]))," ", mean(Array(density_3D[relaxed_indices,:])), std(Array(density_3D[relaxed_indices,:])) )
println("3D temperature slope (t,p,r) = ", mean(Array(temperature_3D)), "+-", std(Array(temperature_3D)),  " ", mean(Array(temperature_3D[perturbed_indices,:])),  "+-", std(Array(temperature_3D[perturbed_indices,:]))," ", mean(Array(temperature_3D[relaxed_indices,:])), std(Array(temperature_3D[relaxed_indices,:])) )
println("3D velocity slope (t,p,r) = ", mean(Array(velocity_3D)), "+-", std(Array(velocity_3D)),  " ", mean(Array(velocity_3D[perturbed_indices,:])),  "+-", std(Array(velocity_3D[perturbed_indices,:]))," ", mean(Array(velocity_3D[relaxed_indices,:])), std(Array(velocity_3D[relaxed_indices,:])) )


println("2D density slope (t,p,r) = ", mean(Array(density)), "+-", std(Array(density)),  " ", mean(Array(density[perturbed_indices,:])),  "+-", std(Array(density[perturbed_indices,:]))," ", mean(Array(density[relaxed_indices,:])), std(Array(density[relaxed_indices,:])) )
println("2D temperature slope (t,p,r) = ", mean(Array(temperature)), "+-", std(Array(temperature)),  " ", mean(Array(temperature[perturbed_indices,:])),  "+-", std(Array(temperature[perturbed_indices,:]))," ", mean(Array(temperature[relaxed_indices,:])), std(Array(temperature[relaxed_indices,:])) )
println("2D velocity slope (t,p,r) = ", mean(Array(velocity)), "+-", std(Array(velocity)),  " ", mean(Array(velocity[perturbed_indices,:])),  "+-", std(Array(velocity[perturbed_indices,:]))," ", mean(Array(velocity[relaxed_indices,:])), std(Array(velocity[relaxed_indices,:])) )


fig, ax = subplots(figsize=(6,6))
ax1 = subplot(2,1,1)
fig.subplots_adjust(top=0.99,bottom=0.125,left=0.165, right=0.99)

x = [minimum(index) - 0.1, maximum(index) + 0.1]
yd = [0.16, 0.16]
yt = [0.02, 0.02]

errorbar(index[perturbed_indices] .- 0.1, density_mean[perturbed_indices], yerr=Array([abs.(density_min[perturbed_indices] - density_mean[perturbed_indices]),
         abs.(density_max[perturbed_indices] - density_mean[perturbed_indices])]), marker=".", markersize=10.0, color="black", linestyle="")
errorbar(index[relaxed_indices] .- 0.1, density_mean[relaxed_indices], yerr=Array([abs.(density_min[relaxed_indices] - density_mean[relaxed_indices]),
         abs.(density_max[relaxed_indices] - density_mean[relaxed_indices])]), marker="d", markersize=6.0, color="black", linestyle="")
plot(x,yd, marker="", color="black", linewidth = 2.0)

errorbar(index[perturbed_indices], temperature_mean[perturbed_indices], yerr=Array([abs.(temperature_min[perturbed_indices] - temperature_mean[perturbed_indices]),
         abs.(temperature_max[perturbed_indices] - temperature_mean[perturbed_indices])]), marker=".", markersize=10.0, color="blue", linestyle="")
errorbar(index[relaxed_indices], temperature_mean[relaxed_indices], yerr=Array([abs.(temperature_min[relaxed_indices] - temperature_mean[relaxed_indices]),
         abs.(temperature_max[relaxed_indices] - temperature_mean[relaxed_indices])]), marker="d", markersize=6.0, color="blue", linestyle="")
plot(x,yt, marker="", color="blue", linewidth = 2.0)

errorbar(index[perturbed_indices] .+ 0.1, velocity_mean[perturbed_indices], yerr=Array([abs.(velocity_min[perturbed_indices] - velocity_mean[perturbed_indices]),
                  abs.(velocity_max[perturbed_indices] - velocity_mean[perturbed_indices])]), marker=".", markersize=10.0, color="red", linestyle="")
errorbar(index[relaxed_indices] .+0.1, velocity_mean[relaxed_indices], yerr=Array([abs.(velocity_min[relaxed_indices] - velocity_mean[relaxed_indices]),
                  abs.(velocity_max[relaxed_indices] - velocity_mean[relaxed_indices])]), marker="d", markersize=6.0, color="red", linestyle="")

xticks(fontsize=11)
yticks(fontsize=11)
ylabel(L"m_{2D}", fontsize = 11)
ax1.set_xticklabels(index,color="w")

ax2 = subplot(2,1,2)

plot(index[perturbed_indices] .- 0.1, density_3D[perturbed_indices], marker=".", markersize=10.0, color="black", linestyle="", label = "Density")
plot(index[relaxed_indices] .- 0.1, density_3D[relaxed_indices], marker="d", markersize=6.0, color="black", linestyle="")



plot(index[perturbed_indices], temperature_3D[perturbed_indices], marker=".", markersize=10.0, color="blue", linestyle="", label = "Temperature")
plot(index[relaxed_indices], temperature_3D[relaxed_indices], marker="d", markersize=6.0, color="blue", linestyle="")


plot(index[perturbed_indices] .+ 0.1, velocity_3D[perturbed_indices], marker=".", markersize=10.0, color="red", linestyle="", label = "Velocity")
plot(index[relaxed_indices] .+ 0.1, velocity_3D[relaxed_indices], marker="d", markersize=6.0, color="red", linestyle="")


xticks(fontsize=11)
yticks(fontsize=11)
ylabel(L"m_{3D}", fontsize = 11)
xlabel("Galaxy cluster", fontsize = 11)
ax2.set_xticklabels(name)
legend(fontsize=11)


plt.show()
