using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit

main = "/home/marco/Scrivania/Ettori_project/Codes/3D_analysis/box_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"

d_t_slope_2d = [0.04,0.09,-0.01,-0.04,-0.03,0.08,-0.06,-0.02]
d_v_slope_2d = [0.67,0.69,0.73,0.54,0.59,0.55,0.68,0.53]
t_v_slope_2d = [0.89,0.80,0.79,0.81,0.62,0.74,0.90,0.65]

d_t_slope_3d = [-0.47,-0.20,-0.45,-0.52,-0.39,-0.57,-0.46,-0.28]
d_v_slope_3d = [1.17,1.37,1.24,1.32,1.27,1.29,1.22,1.24]
t_v_slope_3d = [1.28,1.44,1.33,1.41,1.44,1.37,1.30,1.36]

subplot(3,2,1)
hist(d_t_slope_2d)

xlabel(L"slope(d-t)_{2D}", fontsize=15)
ylabel("N", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)

subplot(3,2,2)
hist(d_t_slope_3d)

xlabel(L"slope(d-t)_{3D}", fontsize=15)
ylabel("N", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)

subplot(3,2,3)
hist(d_v_slope_2d)

xlabel(L"slope(d-v)_{2D}", fontsize=15)
ylabel("N", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)

subplot(3,2,4)
hist(d_v_slope_3d)

xlabel(L"slope(d-v)_{3D}", fontsize=15)
ylabel("N", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)

subplot(3,2,5)
hist(t_v_slope_2d)

xlabel(L"slope(t-v)_{2D}", fontsize=15)
ylabel("N", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)

subplot(3,2,6)
hist(t_v_slope_3d)

xlabel(L"slope(t-v)_{3D}", fontsize=15)
ylabel("N", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)
