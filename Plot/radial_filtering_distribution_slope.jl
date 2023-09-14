using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit


main = "/home/marco/Scrivania/Ettori_project/Codes/3D_analysis/box_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"

d_t_slope_2d = [-0.14,-0.18,-0.07,0.002,-0.43,-0.07,-0.20,-0.12]
d_t_slope_3d = [-0.17,-0.03,-0.05,-0.05,-0.19,-0.01,-0.06,-0.06]

id = [1,2,3,4,5,6,7,8]

mean_2d=(mean(d_t_slope_2d))
mean_3d=(mean(d_t_slope_3d))
x = [1,8]
y_2d = [mean_2d, mean_2d]
y_3d = [mean_3d, mean_3d]
plot(id, d_t_slope_3d, marker=".", markersize="6.0", linestyle="")
plot(x, y_3d)
xlabel("ID", fontsize=15)
ylabel("slope", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)
