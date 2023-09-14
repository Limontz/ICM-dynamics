using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using DataFrames
using LsqFit
using HypothesisTests
@pyimport mpl_toolkits.axes_grid1 as axgrid


main = "/home/marco/Scrivania/Ettori_project/"

file1 = string(main,"Plots/sigmaT-sigman.txt")
file2 = string(main,"Plots/2D_300_box_filtering_mean_fluctuations.txt")

data=readdlm(file1,comments=true)
#sigmaT = data[:,1]
#sigman = data[:,2]
T = data[:,3]

data2=readdlm(file2,comments=true)
sigmaT = data2[:,1]
sigman = data2[:,2]


sigman_perturbed = [data2[8,2],data2[7,2],data2[5,2],data2[4,2],data2[1,2]]
sigmaT_perturbed = [data2[8,1],data2[7,1],data2[5,1],data2[4,1],data2[1,1]]
T_perturbed = [data[8,3],data[7,3],data[5,3],data[4,3],data[1,3]]

sigman_relaxed = [data2[6,2],data2[3,2],data2[2,2]]
sigmaT_relaxed = [data2[6,1],data2[3,1],data2[2,1]]
T_relaxed = [data[6,3],data[3,3],data[2,3]]

scatter(T_perturbed, sigmaT_perturbed, s=30, marker="s", facecolors="blue", edgecolors="black", label="Disturbed")
scatter(T_relaxed, sigmaT_relaxed, s=30, marker="s", facecolors="orange", edgecolors="black", label="Relaxed")
ylabel(L"(\delta_T)_{2D} ", fontsize=20)
xlabel(L"<T>_{2D} [Kev]", fontsize=20)
xticks(fontsize=18)
yticks(fontsize=18)
legend(fontsize=12)
show()
readline()
clf()

scatter(T_perturbed, sigman_perturbed, s=150, marker=".", facecolors="blue", edgecolors="black", label="Disturbed")
scatter(T_relaxed, sigman_relaxed, s=150, marker=".", facecolors="orange", edgecolors="black", label="Relaxed")
ylabel(L"(\delta_{\rho})_{2D}", fontsize=20)
xlabel(L"<T>_{2D} [Kev]", fontsize=20)
xticks(fontsize=18)
yticks(fontsize=18)
legend(fontsize=12)
