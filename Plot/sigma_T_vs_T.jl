using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using DataFrames
using LsqFit

main = "/home/marco/Scrivania/Ettori_project/"

file1=string(main,"IT90_0/01_radial_filtering/IT90_0_2D_temperature_radial_profile.txt")
file2=string(main,"IT90_1/01_radial_filtering/IT90_1_2D_temperature_radial_profile.txt")
file3=string(main,"IT90_2/01_radial_filtering/IT90_2_2D_temperature_radial_profile.txt")
file4=string(main,"IT90_3/01_radial_filtering/IT90_3_2D_temperature_radial_profile.txt")
file5=string(main,"IT90_4/01_radial_filtering/IT90_4_2D_temperature_radial_profile.txt")
file6=string(main,"IT92_0/01_radial_filtering/IT92_0_2D_temperature_radial_profile.txt")
file7=string(main,"IT92_1/01_radial_filtering/IT92_1_2D_temperature_radial_profile.txt")
file8=string(main,"IT92_2/01_radial_filtering/IT92_2_2D_temperature_radial_profile.txt")

T1=readdlm(file1,comments=true)
T2=readdlm(file2,comments=true)
T3=readdlm(file3,comments=true)
T4=readdlm(file4,comments=true)
T5=readdlm(file5,comments=true)
T6=readdlm(file6,comments=true)
T7=readdlm(file7,comments=true)
T8=readdlm(file8,comments=true)

T = cat(T8[:,2],cat(T7[:,2],cat(T6[:,2],cat(T5[:,2],cat(T4[:,2],cat(T3[:,2],cat(T2[:,2],T1[:,2],dims=(1,1)),dims=(1,1)),
              dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))

file1=string(main,"IT90_0/01_radial_filtering/new_IT90_0_2D_radial_filtering_fluctuations_radial_trend.txt")
file2=string(main,"IT90_1/01_radial_filtering/new_IT90_1_2D_radial_filtering_fluctuations_radial_trend.txt")
file3=string(main,"IT90_2/01_radial_filtering/new_IT90_2_2D_radial_filtering_fluctuations_radial_trend.txt")
file4=string(main,"IT90_3/01_radial_filtering/new_IT90_3_2D_radial_filtering_fluctuations_radial_trend.txt")
file5=string(main,"IT90_4/01_radial_filtering/new_IT90_4_2D_radial_filtering_fluctuations_radial_trend.txt")
file6=string(main,"IT92_0/01_radial_filtering/new_IT92_0_2D_radial_filtering_fluctuations_radial_trend.txt")
file7=string(main,"IT92_1/01_radial_filtering/new_IT92_1_2D_radial_filtering_fluctuations_radial_trend.txt")
file8=string(main,"IT92_2/01_radial_filtering/new_IT92_2_2D_radial_filtering_fluctuations_radial_trend.txt")

data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
data3=readdlm(file3,comments=true)
data4=readdlm(file4,comments=true)
data5=readdlm(file5,comments=true)
data6=readdlm(file6,comments=true)
data7=readdlm(file7,comments=true)
data8=readdlm(file8,comments=true)

sigmaT = cat(data8[:,3],cat(data7[:,3],cat(data6[:,3],cat(data5[:,3],cat(data4[:,3],cat(data3[:,3],cat(data2[:,3],data1[:,3],dims=(1,1)),dims=(1,1)),
              dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))

print(T)
plot(T, sigmaT, marker=".", markersize=8.0, linestyle="", color="black")
plt.xlabel(L"T(r)_{2D} / [keV]", fontsize =11)
plt.ylabel(L"\sigma T/T(r)_{2D}", fontsize =11)
plt.yticks(fontsize =11)
plt.xticks(fontsize =11)
plt.show()
