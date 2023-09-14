using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit

d_t_slope_2d_box = [0.03, 0.09, -0.07, -0.07, 0.02, 0.08, -0.25, -0.02]
d_v_slope_2d_box = [0.67, 0.68, 0.78, 0.56, 0.60, 0.56, 0.71, 0.55]
t_v_slope_2d_box = [0.07, 0.08, 0.10, 0.12, 0.06, 0.02, 0.08, 0.05]

d_t_slope_2d_box_error = [0.01, 0.01, 0.01, 0.01, 0.02, 0.01, 0.01, 0.01]
d_v_slope_2d_box_error = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
t_v_slope_2d_box_error = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]

d_t_slope_3d_box = [-0.34, -0.15, -0.35, -0.33, -0.24, -0.42, -0.37, -0.23]
d_v_slope_3d_box = [1.07, 1.25, 1.16, 1.15, 1.14, 1.18, 1.12, 1.09]
t_v_slope_3d_box = [0.45, 0.47, 0.39, 0.46, 0.40, 0.29, 0.35, 0.35]

d_t_slope_3d_box_error = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
d_v_slope_3d_box_error = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
t_v_slope_3d_box_error = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]



d_t_slope_2d_radial = [-0.14, -0.08, -0.07, 0.01, -0.35, -0.11, -0.22, -0.21]
d_v_slope_2d_radial = [1.13, 0.46, 0.61, 0.66, 0.87, 0.80, 0.78, 1.27]
t_v_slope_2d_radial = [0.64, 2.00, 1.08, 1.02, 0.51, 0.45, 1.30, 0.72]

d_t_slope_2d_radial_error = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
d_v_slope_2d_radial_error = [0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
t_v_slope_2d_radial_error = [0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01]

d_t_slope_3d_radial = [-0.16, -0.01, -0.05, -0.05, -0.17, -0.002, -0.07, -0.08]
d_v_slope_3d_radial = [0.99, 0.28, 0.39, 0.82, 0.71, 0.96, 0.27, 0.94]
t_v_slope_3d_radial = [1.33, 2.05, 1.53, 2.61, 1.47, 1.11, 1.48, 1.39]

d_t_slope_3d_radial_error = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
d_v_slope_3d_radial_error = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
t_v_slope_3d_radial_error = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]

#errorbar(d_t_slope_2d_box, d_t_slope_3d_box, xerr=d_t_slope_2d_box_error, yerr=d_t_slope_3d_box_error, marker=".", markersize=10.0, color="darkblue", linestyle="", label="box filtering")
#errorbar(d_t_slope_2d_radial, d_t_slope_3d_radial, xerr=d_t_slope_2d_radial_error, yerr=d_t_slope_3d_radial_error, marker=".", markersize=10.0, color="red", linestyle="", label="radial filtering")
scatter(d_t_slope_2d_box, d_t_slope_3d_box, marker=".", s=250, edgecolor="black", color="darkblue", label="box filtering")
scatter(d_t_slope_2d_radial, d_t_slope_3d_radial, marker=".", s=250, edgecolor="black", color="red", label="radial filtering")
xlabel(L"(\delta T/\delta \rho)_{2D}",fontsize=20)
ylabel(L"(\delta T/\delta \rho)_{3D}",fontsize=20)
xlim(-0.5, 0.2)
ylim(-0.5, 0.2)
xticks(fontsize=18)
yticks(fontsize=18)
legend(fontsize=15, loc="upper left")
show()

readline()
clf()

#errorbar(d_v_slope_2d_box, d_v_slope_3d_box, xerr=d_v_slope_2d_box_error, yerr=d_v_slope_3d_box_error, marker=".", markersize=10.0, color="darkblue", linestyle="", label="box filtering")
#errorbar(d_v_slope_2d_radial, d_v_slope_3d_radial, xerr=d_v_slope_2d_radial_error, yerr=d_v_slope_3d_radial_error, marker=".", markersize=10.0, color="red", linestyle="", label="radial filtering")
scatter(d_v_slope_2d_box, d_v_slope_3d_box, marker="s", s=50, edgecolor="black", color="darkblue")
scatter(d_v_slope_2d_radial, d_v_slope_3d_radial, marker="s", s=50, edgecolor="black", color="red")
xlabel(L"(\delta v/\delta \rho)_{2D}",fontsize=20)
ylabel(L"(\delta v/\delta \rho)_{3D}",fontsize=20)
xlim(0.2, 1.3)
ylim(0.2, 1.3)
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
xticks(fontsize=18)
yticks(fontsize=18)
#legend(fontsize=15, loc="upper left")
show()
