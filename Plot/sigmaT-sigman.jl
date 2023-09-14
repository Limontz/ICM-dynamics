using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using DataFrames
using LsqFit
using HypothesisTests
@pyimport mpl_toolkits.axes_grid1 as axgrid


main = "/home/marco/Scrivania/Ettori_project/"

file1 = string(main,"Plots/x_sigmaT-sigman_2D.txt")
file2 = string(main,"Plots/y_sigmaT-sigman_2D.txt")
file3 = string(main,"Plots/z_sigmaT-sigman_2D.txt")
file4 = string(main,"Plots/sigmaT-sigman.txt")

data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
data3=readdlm(file3,comments=true)
data4 = readdlm(file4,comments=true)

df = DataFrame( x = data1[:,1]./data1[:,3],
                y = data2[:,1]./data2[:,3],
                z = data3[:,1]./data3[:,3] )

T = DataFrame( Tx = data1[:,3], Ty = data2[:,3], Tz = data3[:,3])


fluct_mean = (mean(Array(df[:,:]), dims=2))
fluct_max, fluct_min = (maximum(Array(df[:,:]), dims=2)), (minimum(Array(df[:,:]), dims=2))

T_mean = (mean(Array(T[:,:]), dims=2))
T_max, T_min = (maximum(Array(T[:,:]), dims=2)), (minimum(Array(T[:,:]), dims=2))

fluct_3D = data4[:,1]./data4[:,3]
T_3D = data4[:,3]

perturbed_indices = [1,4,5,7,8]
relaxed_indices = [2,3,6]

#Lovisari file

file1 = string(main,"/scatter_R500_Tne_sn30.txt")
lovisaridata=readdlm(file1,comments=true)

fig, ax = subplots(figsize=(6,6))
fig.subplots_adjust(top=0.99,bottom=0.125,left=0.165, right=0.99)

errorbar(T_mean[perturbed_indices], fluct_mean[perturbed_indices],
         yerr=Array([abs.(fluct_min[perturbed_indices] - fluct_mean[perturbed_indices]),
         abs.(fluct_max[perturbed_indices] - fluct_mean[perturbed_indices])]),
         xerr=Array([abs.(T_min[perturbed_indices] - T_mean[perturbed_indices]),
         abs.(T_max[perturbed_indices] - T_mean[perturbed_indices])]), linestyle="",
         marker=".", markersize="8.0", color="blue", label="Disturbed")

errorbar(T_mean[relaxed_indices], fluct_mean[relaxed_indices],
         yerr=Array([abs.(fluct_min[relaxed_indices] - fluct_mean[relaxed_indices]),
         abs.(fluct_max[relaxed_indices] - fluct_mean[relaxed_indices])]),
         xerr=Array([abs.(T_min[relaxed_indices] - T_mean[relaxed_indices]),
         abs.(T_max[relaxed_indices] - T_mean[relaxed_indices])]), linestyle="",
         marker=".", markersize="8.0", color="orange", label="Relaxed")

plot(lovisaridata[:,4], lovisaridata[:,6]./lovisaridata[:,4], marker = ".", color="darkgreen", linestyle="",
     markersize="8.0", label="Lovisari(2023)")

ylabel(L"(\sigma_{T}/<T>)_{2D} ", fontsize=16)
xlabel(L"<T>_{2D} [Kev]", fontsize=16)
xticks([0.6, 1 ,2, 5, 10],fontsize=13)
yticks(fontsize=13)
legend(fontsize=13)
xscale("log")
xticks([0.6, 1, 2, 5, 10], ["0.6", "1", "2", "5", "10"], fontsize=13)
show()
readline()
clf()

fig, ax = subplots(figsize=(6,6))
fig.subplots_adjust(top=0.99,bottom=0.125,left=0.165, right=0.99)

plot(T_3D[perturbed_indices], fluct_3D[perturbed_indices], marker=".", markersize = "8.0", color="blue", label="Disturbed", linestyle="")
plot(T_3D[relaxed_indices], fluct_3D[relaxed_indices], marker=".", markersize="8.0",color="orange", label="Relaxed", linestyle="")
ylabel(L"(\sigma_{T}/<T>)_{3D}", fontsize=16)
xlabel(L"<T>_{3D} [Kev]", fontsize=16)
xticks(fontsize=13)
yticks([0.30, 0.35, 0.40, 0.45, 0.50],fontsize=13)
