using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using DataFrames


main = "/home/marco/Scrivania/Ettori_project/"

file1=string(main,"IT90_0/01_radial_filtering/new_IT90_0_3D_radial_filtering_slope_radial_trend.txt")
file2=string(main,"IT90_1/01_radial_filtering/new_IT90_1_3D_radial_filtering_slope_radial_trend.txt")
file3=string(main,"IT90_2/01_radial_filtering/new_IT90_2_3D_radial_filtering_slope_radial_trend.txt")
file4=string(main,"IT90_3/01_radial_filtering/new_IT90_3_3D_radial_filtering_slope_radial_trend.txt")
file5=string(main,"IT90_4/01_radial_filtering/new_IT90_4_3D_radial_filtering_slope_radial_trend.txt")
file6=string(main,"IT92_0/01_radial_filtering/new_IT92_0_3D_radial_filtering_slope_radial_trend.txt")
file7=string(main,"IT92_1/01_radial_filtering/new_IT92_1_3D_radial_filtering_slope_radial_trend.txt")
file8=string(main,"IT92_2/01_radial_filtering/new_IT92_2_3D_radial_filtering_slope_radial_trend.txt")

data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
data3=readdlm(file3,comments=true)
data4=readdlm(file4,comments=true)
data5=readdlm(file5,comments=true)
data6=readdlm(file6,comments=true)
data7=readdlm(file7,comments=true)
data8=readdlm(file8,comments=true)

slope_d_t = DataFrame(bin = data2[:,1], it900=data1[:,2],
                    it901 = data2[:,2],it902=data3[:,2],
                    it903=data4[:,2],it904=data5[:,2],
                    it920=data6[:,2],it921=data7[:,2],
                    it922=data8[:,2])

slope_d_v = DataFrame(bin = data2[:,1], it900=data1[:,3],
                    it901 = data2[:,3],it902=data3[:,3],
                    it903=data4[:,3],it904=data5[:,3],
                    it920=data6[:,3],it921=data7[:,3],
                    it922=data8[:,3])

slope_t_v = DataFrame(bin = data2[:,1], it900=data1[:,5],
                    it901 = data2[:,5],it902=data3[:,5],
                    it903=data4[:,5],it904=data5[:,5],
                    it920=data6[:,5],it921=data7[:,5],
                    it922=data8[:,5])

println(slope_d_t)
d_t_mean, d_t_sigma=(mean(Array(slope_d_t[:,2:9]), dims=2)), std(Array(slope_d_t[:,2:9]), dims=2)


d_v_mean, d_v_sigma=(mean(Array(slope_d_v[:,2:9]), dims=2)), std(Array(slope_d_v[:,2:9]), dims=2)
t_v_mean, t_v_sigma=(mean(Array(slope_t_v[:,2:9]), dims=2)), std(Array(slope_t_v[:,2:9]), dims=2)

fig, ax1 = subplots(figsize=(6,8))
subplot(3,1,1)
fig.subplots_adjust(top=0.99,bottom=0.125,left=0.165, right=0.99)
for i in 2:length(slope_d_t[1,:])
    plot(slope_d_t.bin/10, slope_d_t[:,i], marker=".", markersize=10.0, linestyle="", color="black")
end

errorbar(slope_d_t.bin/10, d_t_mean, yerr=d_t_sigma, marker=".", markersize=10.0, color="red", linestyle="")
frame = gca()
frame.axes.xaxis.set_ticklabels([])
ylabel(L"(\delta T / \delta \rho)_{3D}",fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)

subplot(3,1,2)
for i in 2:length(slope_d_v[1,:])
    plot(slope_d_v.bin/10, slope_d_v[:,i], marker="s", markersize=6.0, linestyle="", color="black")
end
errorbar(slope_d_v.bin/10, d_v_mean, yerr=d_v_sigma, marker="s", markersize=6.0, color="red", linestyle="")
ylabel(L"(\delta v / \delta \rho)_{3D}",fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)
frame = gca()
frame.axes.xaxis.set_ticklabels([])


subplot(3,1,3)
for i in 2:length(slope_t_v[1,:])
    plot(slope_t_v.bin/10, slope_t_v[:,i], marker="D", markersize=5.0, linestyle="", color="black")
end
errorbar(slope_t_v.bin/10, t_v_mean, yerr=t_v_sigma, marker="D", markersize=5.0, color="red", linestyle="")
ylabel(L"(\delta v^2 / \delta T)_{3D}",fontsize=15)
xlabel(L"R/R_{500}",fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)
show()
