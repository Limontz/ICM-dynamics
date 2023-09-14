using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit


main = "/home/marco/Scrivania/Ettori_project/Codes/3D_analysis/box_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"

d_t_slope_2d = [-0.14,-0.18,-0.07,0.002,-0.43,-0.07,-0.20,-0.12]
d_v_slope_2d = [1.13,0.78,0.62,0.67,0.92,0.87,0.77,1.25]
t_v_slope_2d = [1.52,1.83,1.48,1.89,1.06,1.60,1.85,1.48]

d_t_slope_3d = [-0.17,-0.03,-0.05,-0.05,-0.19,-0.01,-0.06,-0.06]
d_v_slope_3d = [1.0,0.36,0.39,0.83,0.72,1.03,0.28,0.92]
t_v_slope_3d = [1.53,1.91,1.65,2.24,1.63,1.73,1.79,1.62]

plot(d_t_slope_3d, d_t_slope_2d, marker=".", markersize="6.0", linestyle="", color="black", label="dens-temp")
#plot(d_v_slope_3d, d_v_slope_2d, marker=".", markersize="6.0", linestyle="", color="blue", label="dens-vel")
#plot(t_v_slope_3d, t_v_slope_2d, marker=".", markersize="6.0", linestyle="", color="red", label="temp-vel")
xlabel(L"3D slope",fontsize=15)
ylabel(L"2D slope", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)
legend()
