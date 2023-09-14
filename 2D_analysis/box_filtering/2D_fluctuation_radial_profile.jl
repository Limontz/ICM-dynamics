using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using DataFrames


main = "/home/marco/Scrivania/Ettori_project/"

#file1=string(main,"IT90_0/IT90_0_2D_box_filtering_fluctuations_radial_trend.txt")
file2=string(main,"IT90_1/IT90_1_2D_box_filtering_fluctuations_radial_trend.txt")
#file3=string(main,"IT90_2/IT90_2_2D_box_filtering_fluctuations_radial_trend.txt")
#file4=string(main,"IT90_3/IT90_3_2D_box_filtering_fluctuations_radial_trend.txt")
#file5=string(main,"IT90_4/IT90_4_2D_box_filtering_fluctuations_radial_trend.txt")
#file6=string(main,"IT92_0/IT92_0_2D_box_filtering_fluctuations_radial_trend.txt")
#file7=string(main,"IT92_1/IT92_1_2D_box_filtering_fluctuations_radial_trend.txt")
#file8=string(main,"IT92_2/IT92_2_2D_box_filtering_fluctuations_radial_trend.txt")

#data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
#data3=readdlm(file3,comments=true)
#data4=readdlm(file4,comments=true)
#data5=readdlm(file5,comments=true)
#data6=readdlm(file6,comments=true)
#data7=readdlm(file7,comments=true)
#data8=readdlm(file8,comments=true)

density = DataFrame(bin = data2[:,1],
                    it901 = data2[:,2])

temperature = DataFrame(bin = data2[:,1],
                    it901 = data2[:,3])

velocity = DataFrame(bin = data2[:,1],
                    it901 = data2[:,4])

print(density.bin)
print(density.it901)
density_mean, density_sigma=(mean(Array(density[:,2]), dims=2)), std(Array(density[:,1:2]), dims=2)
temperature_mean, temperature_sigma=(mean(Array(temperature[:,2]), dims=2)), std(Array(temperature[:,1:2]), dims=2)
velocity_mean, velocity_sigma=(mean(Array(velocity[:,2]), dims=2)), std(Array(velocity[:,1:2]), dims=2)

for i in 2:length(density[1,:])
    plot(density.bin, density[:,i], marker=".", markersize=8.0, linestyle="", color="black")
end
errorbar(density.bin, density_mean, yerr=density_sigma, marker=".", markersize=8.0, color="red", linestyle="")
show()
