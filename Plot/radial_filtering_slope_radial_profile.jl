using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using DataFrames
using LsqFit

main = "/home/marco/Scrivania/Ettori_project/"

file1=string(main,"IT90_0/02_radial_filtering/IT90_0_2D_radial_filtering_slope_radial_trend.txt")
file2=string(main,"IT90_1/02_radial_filtering/IT90_1_2D_radial_filtering_slope_radial_trend.txt")
file3=string(main,"IT90_2/02_radial_filtering/IT90_2_2D_radial_filtering_slope_radial_trend.txt")
file4=string(main,"IT90_3/02_radial_filtering/IT90_3_2D_radial_filtering_slope_radial_trend.txt")
file5=string(main,"IT90_4/02_radial_filtering/IT90_4_2D_radial_filtering_slope_radial_trend.txt")
file6=string(main,"IT92_0/02_radial_filtering/IT92_0_2D_radial_filtering_slope_radial_trend.txt")
file7=string(main,"IT92_1/02_radial_filtering/IT92_1_2D_radial_filtering_slope_radial_trend.txt")
file8=string(main,"IT92_2/02_radial_filtering/IT92_2_2D_radial_filtering_slope_radial_trend.txt")


data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
data3=readdlm(file3,comments=true)
data4=readdlm(file4,comments=true)
data5=readdlm(file5,comments=true)
data6=readdlm(file6,comments=true)
data7=readdlm(file7,comments=true)
data8=readdlm(file8,comments=true)

density_2D = DataFrame(bin = data2[:,1], it900=data1[:,2],
                    it901 = data2[:,2],it902=data3[:,2],
                    it903=data4[:,2],it904=data5[:,2],
                    it920=data6[:,2],it921=data7[:,2],
                    it922=data8[:,2])

temperature_2D = DataFrame(bin = data2[:,1], it900=data1[:,3],
                    it901 = data2[:,3],it902=data3[:,3],
                    it903=data4[:,3],it904=data5[:,3],
                    it920=data6[:,3],it921=data7[:,3],
                    it922=data8[:,3])

velocity_2D = DataFrame(bin = data2[:,1], it900=data1[:,4],
                    it901 = data2[:,4],it902=data3[:,4],
                    it903=data4[:,4],it904=data5[:,4],
                    it920=data6[:,4],it921=data7[:,4],
                    it922=data8[:,4])

density_mean_2D, density_sigma_2D=(mean(Array(density_2D[:,2:9]), dims=2)), std(Array(density_2D[:,2:9]), dims=2)
temperature_mean_2D, temperature_sigma_2D=(mean(Array(temperature_2D[:,2:9]), dims=2)), std(Array(temperature_2D[:,2:9]), dims=2)
velocity_mean_2D, velocity_sigma_2D=(mean(Array(velocity_2D[:,2:9]), dims=2)), std(Array(velocity_2D[:,2:9]), dims=2)


file1=string(main,"IT90_0/02_radial_filtering/new_IT90_0_3D_radial_filtering_slope_radial_trend.txt")
file2=string(main,"IT90_1/02_radial_filtering/new_IT90_1_3D_radial_filtering_slope_radial_trend.txt")
file3=string(main,"IT90_2/02_radial_filtering/new_IT90_2_3D_radial_filtering_slope_radial_trend.txt")
file4=string(main,"IT90_3/02_radial_filtering/new_IT90_3_3D_radial_filtering_slope_radial_trend.txt")
file5=string(main,"IT90_4/02_radial_filtering/new_IT90_4_3D_radial_filtering_slope_radial_trend.txt")
file6=string(main,"IT92_0/02_radial_filtering/new_IT92_0_3D_radial_filtering_slope_radial_trend.txt")
file7=string(main,"IT92_1/02_radial_filtering/new_IT92_1_3D_radial_filtering_slope_radial_trend.txt")
file8=string(main,"IT92_2/02_radial_filtering/new_IT92_2_3D_radial_filtering_slope_radial_trend.txt")


data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
data3=readdlm(file3,comments=true)
data4=readdlm(file4,comments=true)
data5=readdlm(file5,comments=true)
data6=readdlm(file6,comments=true)
data7=readdlm(file7,comments=true)
data8=readdlm(file8,comments=true)

density_3D = DataFrame(bin = data2[:,1], it900=data1[:,2],
                    it901 = data2[:,2],it902=data3[:,2],
                    it903=data4[:,2],it904=data5[:,2],
                    it920=data6[:,2],it921=data7[:,2],
                    it922=data8[:,2])

temperature_3D = DataFrame(bin = data2[:,1], it900=data1[:,3],
                    it901 = data2[:,3],it902=data3[:,3],
                    it903=data4[:,3],it904=data5[:,3],
                    it920=data6[:,3],it921=data7[:,3],
                    it922=data8[:,3])

velocity_3D = DataFrame(bin = data2[:,1], it900=data1[:,4],
                    it901 = data2[:,4],it902=data3[:,4],
                    it903=data4[:,4],it904=data5[:,4],
                    it920=data6[:,4],it921=data7[:,4],
                    it922=data8[:,4])


density_mean_3D, density_sigma_3D=(mean(Array(density_3D[:,2:9]), dims=2)), std(Array(density_3D[:,3:9]), dims=2)
temperature_mean_3D, temperature_sigma_3D=(mean(Array(temperature_3D[:,2:9]), dims=2)), std(Array(temperature_3D[:,2:9]), dims=2)
velocity_mean_3D, velocity_sigma_3D=(mean(Array(velocity_3D[:,2:9]), dims=2)), std(Array(velocity_3D[:,2:9]), dims=2)


fig, ax1 = subplots(figsize=(6,8))
subplot(3,1,1)
fig.subplots_adjust(top=0.99,bottom=0.125,left=0.165, right=0.99)

errorbar(density_2D.bin/10, density_mean_2D, yerr=density_sigma_2D, marker=".", markersize=8.0, color="blue", linestyle="", label="2D")
errorbar(density_2D.bin/10 .+ 0.05, density_mean_3D, yerr=density_sigma_3D, marker=".", markersize=8.0, color="red", linestyle="", label="3D")
ylabel(L"\delta T/ \delta \rho",fontsize=13)
frame = gca()
frame.axes.xaxis.set_ticklabels([])
xticks(fontsize=11)
yticks(fontsize=11)
plt.legend()

subplot(3,1,2)

errorbar(temperature_2D.bin/10, temperature_mean_2D, yerr=temperature_sigma_2D, marker="s", markersize=5.0, color="blue", linestyle="")
errorbar(temperature_3D.bin/10 .+ 0.05, temperature_mean_3D, yerr=temperature_sigma_3D, marker="s", markersize=5.0, color="red", linestyle="")

ylabel(L"\delta v/ \delta \rho",fontsize=13)
xticks(fontsize=11)
yticks(fontsize=11)
frame = gca()
frame.axes.xaxis.set_ticklabels([])


subplot(3,1,3)

errorbar(velocity_2D.bin/10, velocity_mean_2D, yerr=velocity_sigma_2D, marker="D", markersize=4.0, color="blue", linestyle="")
errorbar(velocity_3D.bin/10 .+ 0.05, velocity_mean_3D, yerr=velocity_sigma_3D, marker="D", markersize=4.0, color="red", linestyle="")
xlabel(L"R/R_{500}",fontsize=11)
ylabel(L"\delta v/ \delta T",fontsize=13)
xticks(fontsize=11)
yticks(fontsize=11)
show()
