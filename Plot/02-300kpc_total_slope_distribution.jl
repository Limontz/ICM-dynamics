using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit

d_t_slope_2d_box = [0.02, -0.01, -0.13, -0.11, -0.06, -0.05, -0.28, -0.09]
d_v_slope_2d_box = [0.67, 0.75, 0.77, 0.62, 0.64, 0.72, 1.10, 0.57]
t_v_slope_2d_box = [0.21, 0.22, 0.39, 0.55, 0.21, 0.08, 0.84, 0.16]

#d_t_slope_2d_box_error = []
#d_v_slope_2d_box_error = []
#t_v_slope_2d_box_error = []

d_t_slope_3d_box = [-0.30, -0.04, -0.29, -0.25, -0.23, -0.36, -0.29, -0.12]
d_v_slope_3d_box = [1.25, 1.53, 1.25, 1.28, 1.29, 1.45, 0.69, 1.24]
t_v_slope_3d_box = [1.33, 1.12, 0.92, 1.52, 1.06, 0.73, 0.26, 0.98]

#d_t_slope_3d_box_error = []
#d_v_slope_3d_box_error = []
#t_v_slope_3d_box_error = []



d_t_slope_2d_radial = [0.05, -0.06, -0.04, 0.05, -0.18, 0.05, -0.16, -0.01]
d_v_slope_2d_radial = [0.82, 0.44, 0.55, 0.57, 0.73, 0.59, 0.66, 0.69]
t_v_slope_2d_radial = [0.55, 1.78, 0.95, 0.86, 0.46,0.35, 1.19, 0.59]

#d_t_slope_2d_radial_error = []
#d_v_slope_2d_radial_error = []
#t_v_slope_2d_radial_error = []

d_t_slope_3d_radial = [-0.11, -0.01, -0.05, -0.02, -0.14, 0.01, -0.06, -0.05]
d_v_slope_3d_radial = [0.88, 0.27, 0.39, 0.75, 0.67, 0.92, 0.27, 0.86]
t_v_slope_3d_radial = [1.32, 2.01, 1.52, 2.52, 1.45, 1.11, 1.47, 1.38]

#d_t_slope_3d_radial_error = []
#d_v_slope_3d_radial_error = []
#t_v_slope_3d_radial_error = []


scatter(d_t_slope_2d_box, d_t_slope_3d_box, marker=".", s=250, edgecolor="black", color="darkblue", label="box filtering")
scatter(d_t_slope_2d_radial, d_t_slope_3d_radial, marker=".", s=250, edgecolor="black", color="red", label="radial filtering")
xlabel(L"(\delta T/\delta \rho)_{2D}",fontsize=20)
ylabel(L"(\delta T/\delta \rho)_{3D}",fontsize=20)
#xlim(-0.4, 0.2)
#ylim(-0.4, 0.2)
xticks([-0.30, -0.20,-0.10,0], fontsize=18)
yticks([-0.30, -0.20,-0.10,0], fontsize=18)
legend(fontsize=15, loc="lower left")
show()

readline()
clf()

#errorbar(d_v_slope_2d_box, d_v_slope_3d_box, xerr=d_v_slope_2d_box_error, yerr=d_v_slope_3d_box_error, marker=".", markersize=10.0, color="darkblue", linestyle="", label="box filtering")
#errorbar(d_v_slope_2d_radial, d_v_slope_3d_radial, xerr=d_v_slope_2d_radial_error, yerr=d_v_slope_3d_radial_error, marker=".", markersize=10.0, color="red", linestyle="", label="radial filtering")
scatter(d_v_slope_2d_box, d_v_slope_3d_box, marker="s", s=50, edgecolor="black", color="darkblue")
scatter(d_v_slope_2d_radial, d_v_slope_3d_radial, marker="s", s=50, edgecolor="black", color="red")
xlabel(L"(\delta v/\delta \rho)_{2D}",fontsize=20)
ylabel(L"(\delta v/\delta \rho)_{3D}",fontsize=20)
xlim(0.2, 1.2)
ylim(0.2, 1.7)
xticks(fontsize=18)
yticks(fontsize=18)
#legend(fontsize=15, loc="upper left")
show()

readline()
clf()


scatter(t_v_slope_2d_box, t_v_slope_3d_box, marker="D", s=50, edgecolor="black", color="darkblue")
scatter(t_v_slope_2d_radial, t_v_slope_3d_radial, marker="D", s=50, edgecolor="black", color="red")
xlabel(L"(\delta v/\delta T)_{2D}",fontsize=20)
ylabel(L"(\delta v/\delta T)_{3D}",fontsize=20)
#xlim(0.2, 1.3)
#ylim(0.2, 1.3)
xticks([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8],fontsize=18)
yticks(fontsize=18)
#legend(fontsize=15, loc="upper left")
show()
