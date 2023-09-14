using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit

main = "/home/marco/Scrivania/Tesi/Ammassi_tesi/"
snap = "0242"

const mu = 0.62
const mp = 1.67e-24
const kb = 1.28e-16
const pc = 3.08e18
const kpc = (10.0^3) * pc
const G = 6.67e-8
const gamma = 5.0 / 3.0

rvir=1.01*1.e3*kpc
r500 = 5*rvir/7

include(string(main, "read.jl"))
d, d_dm, temp, vx, vy, vz, Npoint, conv_d  = read(snap, main)
Nhalf = 150

d_max, index = findmax(d)
xc = index[1]
yc = index[2]
zc = index[3]

#*********************PROJECTION***********************************

include(string(main, "Projection.jl"))

d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D = projection(d, temp, vx, vy, vz, Npoint)
#****************** RADIAL PROFILE ****************************
include(string(main, "radial_profile.jl"))

r_max = 20. * kpc * Nhalf * sqrt(3)
rmin = 20. * kpc
radial_bin = 300
temp_prof, radius =profile(temp, xc, yc, zc, r_max, rmin, Nhalf,radial_bin)##

include(string(main, "2D_radial_profile.jl"))
radial_bin = 320
d_2D_prof, temp_2D_prof, radius = profile2D(d_2D, spec_T_2D, xc, yc, zc, radial_bin)


#****************** FILTERING *********************************

include(string(main, "v_turb2.jl"))
include(string(main, "v_mean2.jl"))
deltav, deltad, deltaT, L = v_turb2(d, vx, vy, vz, temp, temp_prof, Npoint)

include(string(main, "2D_filtering.jl"))
include(string(main, "v_mean2D.jl"))
deltad_2D, deltaT_2D, deltav_2D = filtering2D(d_2D, d_2D_prof, spec_T_2D, temp_2D_prof,
                                   vx_2D, vy_2D, vz_2D, radial_bin, Npoint)


dx=0.02
Ncell = floor(Int32, (r500/(20*kpc)))

ratio_3D = zeros(2*Ncell,2*Ncell)
ratio_2D = zeros(2*Ncell,2*Ncell)
@inbounds for i in 1:2*Ncell
@inbounds @simd for j in 1:2*Ncell
                    ratio_3D[i,j] = deltaT[xc-Ncell+i,yc-Ncell+j,zc]/deltad[xc-Ncell+i,yc-Ncell+j,zc]
                    if (ratio_3D[i,j] > 1)
                        ratio_3D[i,j] = 1
                    end
                    if (ratio_3D[i,j] < -1)
                        ratio_3D[i,j] = -1
                    end
                    ratio_2D[i,j] = deltaT_2D[xc-Ncell+i,yc-Ncell+j]/deltad_2D[xc-Ncell+i,yc-Ncell+j]

                    if (ratio_2D[i,j] > 1)
                        ratio_2D[i,j] = 1
                    end
                    if (ratio_2D[i,j] < -1)
                        ratio_2D[i,j] = -1
                    end
                end
          end

mycmap = matplotlib.colors.ListedColormap(["red", "white", "blue"])
bounds = [-1.,-0.1,0.1,1.]
subplot(3,2,1)
pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((deltaT[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc]))))

axis("scaled")
#title("IT90_3", fontsize= "20")
xticks(fontsize= 15)
yticks(fontsize= 15)
cbar=colorbar(orientation="vertical")
cbar.ax.tick_params(labelsize=15)
cbar.ax.yaxis.get_offset_text().set(size=15)
#cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"\Delta T / T (3D)",fontsize= 15)
xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 15)
ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 15)

subplot(3,2,2)
pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((deltaT_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]))))

axis("scaled")
#title("IT90_3", fontsize= "20")
xticks(fontsize= 15)
yticks(fontsize= 15)
cbar=colorbar(orientation="vertical")
cbar.ax.tick_params(labelsize=15)
cbar.ax.yaxis.get_offset_text().set(size=15)
#cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"\Delta T / T (2D)",fontsize= 15)
xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 15)
ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 15)

subplot(3,2,3)
pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((deltad[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc]))))

axis("scaled")
#title("IT90_3", fontsize= "20")
xticks(fontsize= 15)
yticks(fontsize= 15)
cbar=colorbar(orientation="vertical")
cbar.ax.tick_params(labelsize=15)
cbar.ax.yaxis.get_offset_text().set(size=15)
#cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"\Delta \rho / \rho (3D)",fontsize= 15)
xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 15)
ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 15)

subplot(3,2,4)
pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((deltad_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]))))

axis("scaled")
#title("IT90_3", fontsize= "20")
xticks(fontsize= 15)
yticks(fontsize= 15)
cbar=colorbar(orientation="vertical")
cbar.ax.tick_params(labelsize=15)
cbar.ax.yaxis.get_offset_text().set(size=15)
#cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"\Delta \rho / ρ (2D)",fontsize= 15)
xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 15)
ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 15)



subplot(3,2,5)
pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((ratio_3D))), cmap=mycmap)

axis("scaled")
#title("IT90_3", fontsize= "20")
xticks(fontsize= 15)
yticks(fontsize= 15)
cbar=colorbar(boundaries=bounds,ticks=bounds,orientation="vertical")
cbar.ax.tick_params(labelsize=15)
cbar.ax.yaxis.get_offset_text().set(size=15)
#cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"G-1 (3D)",fontsize= 15)
xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 15)
ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 15)

subplot(3,2,6)
pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((ratio_2D))), cmap=mycmap)

axis("scaled")
#title("IT90_3", fontsize= "20")
xticks(fontsize= 15)
yticks(fontsize= 15)
cbar=colorbar(boundaries=bounds, ticks=bounds, orientation="vertical")
cbar.ax.tick_params(labelsize=15)
cbar.ax.yaxis.get_offset_text().set(size=15)
#cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"G-1 (3D)",fontsize= 15)
xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 15)
ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 15)
