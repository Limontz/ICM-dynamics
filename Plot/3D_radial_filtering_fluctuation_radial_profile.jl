using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using DataFrames


main = "/home/marco/Scrivania/Ettori_project/"

file1=string(main,"IT90_0/01_radial_filtering/IT90_0_3D_radial_filtering_fluctuations_radial_trend.txt")
file2=string(main,"IT90_1/01_radial_filtering/IT90_1_3D_radial_filtering_fluctuations_radial_trend.txt")
file3=string(main,"IT90_2/01_radial_filtering/IT90_2_3D_radial_filtering_fluctuations_radial_trend.txt")
file4=string(main,"IT90_3/01_radial_filtering/IT90_3_3D_radial_filtering_fluctuations_radial_trend.txt")
file5=string(main,"IT90_4/01_radial_filtering/IT90_4_3D_radial_filtering_fluctuations_radial_trend.txt")
file6=string(main,"IT92_0/01_radial_filtering/IT92_0_3D_radial_filtering_fluctuations_radial_trend.txt")
file7=string(main,"IT92_1/01_radial_filtering/IT92_1_3D_radial_filtering_fluctuations_radial_trend.txt")
file8=string(main,"IT92_2/01_radial_filtering/IT92_2_3D_radial_filtering_fluctuations_radial_trend.txt")

data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
data3=readdlm(file3,comments=true)
data4=readdlm(file4,comments=true)
data5=readdlm(file5,comments=true)
data6=readdlm(file6,comments=true)
data7=readdlm(file7,comments=true)
data8=readdlm(file8,comments=true)

density = DataFrame(bin = data2[:,1], it900=data1[:,2],
                    it901 = data2[:,2],it902=data3[:,2],
                    it903=data4[:,2],it904=data5[:,2],
                    it920=data6[:,2],it921=data7[:,2],
                    it922=data8[:,2])

temperature = DataFrame(bin = data2[:,1], it900=data1[:,3],
                    it901 = data2[:,3],it902=data3[:,3],
                    it903=data4[:,3],it904=data5[:,3],
                    it920=data6[:,3],it921=data7[:,3],
                    it922=data8[:,3])

velocity = DataFrame(bin = data2[:,1], it900=data1[:,4],
                    it901 = data2[:,4],it902=data3[:,4],
                    it903=data4[:,4],it904=data5[:,4],
                    it920=data6[:,4],it921=data7[:,4],
                    it922=data8[:,4])

density_mean, density_sigma=(mean(Array(density[:,2:9]), dims=2)), std(Array(density[:,3:9]), dims=2)
println(density)
exit()

temperature_mean, temperature_sigma=(mean(Array(temperature[:,2:9]), dims=2)), std(Array(temperature[:,2:9]), dims=2)
velocity_mean, velocity_sigma=(mean(Array(velocity[:,2:9]), dims=2)), std(Array(velocity[:,2:9]), dims=2)

fig, ax1 = subplots(figsize=(6,8))
subplot(3,1,1)
fig.subplots_adjust(top=0.99,bottom=0.125,left=0.165, right=0.99)
for i in 2:length(density[1,:])
    plot(density.bin/10, density[:,i], marker=".", markersize=8.0, linestyle="", color="black")
end

errorbar(density.bin/10, density_mean, yerr=density_sigma, marker=".", markersize=8.0, color="red", linestyle="")
ylabel(L"\sigma_{\rho,3D}",fontsize=13)
frame = gca()
frame.axes.xaxis.set_ticklabels([])
xticks(fontsize=11)
yticks(fontsize=11)

subplot(3,1,2)
for i in 2:length(temperature[1,:])
    plot(temperature.bin/10, temperature[:,i], marker=".", markersize=8.0, linestyle="", color="black")
end
errorbar(temperature.bin/10, temperature_mean, yerr=temperature_sigma, marker=".", markersize=8.0, color="red", linestyle="")
ylabel(L"\sigma_{T,3D}",fontsize=13)
xticks(fontsize=11)
yticks(fontsize=11)
frame = gca()
frame.axes.xaxis.set_ticklabels([])


subplot(3,1,3)
for i in 2:length(velocity[1,:])
    plot(velocity.bin/10, velocity[:,i], marker=".", markersize=8.0, linestyle="", color="black")
end
errorbar(velocity.bin/10, velocity_mean, yerr=velocity_sigma, marker=".", markersize=8.0, color="red", linestyle="")
xlabel(L"R/R_{500}",fontsize=11)
ylabel(L"\sigma_{v,3D} ",fontsize=13)
xticks(fontsize=11)
yticks(fontsize=11)
show()
