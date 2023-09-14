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
