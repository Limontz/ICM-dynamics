using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit


main = "/home/marco/Scrivania/Ettori_project/Codes/3D_analysis/box_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"

d_t_slope_2d = [0.04,0.09,-0.01,-0.04,-0.03,0.08,-0.06,-0.02]
d_t_slope_3d = [-0.47,-0.20,-0.45,-0.52,-0.39,-0.57,-0.46,-0.28]

id = [1,2,3,4,5,6,7,8]

mean_2d=(mean(d_t_slope_2d))
mean_3d=(mean(d_t_slope_3d))
x = [1,8]
y_2d = [mean_2d, mean_2d]
y_3d = [mean_3d, mean_3d]
plot(id, d_t_slope_2d, marker=".", markersize="6.0", linestyle="")
plot(x, y_2d)
xlabel("ID", fontsize=15)
ylabel("slope", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)
