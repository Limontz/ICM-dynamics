using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit


main = "/home/marco/Scrivania/Ettori_project/Codes/3D_analysis/box_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"
snap = "0187"

const mu = 0.62
const mp = 1.67e-24
const kb = 1.28e-16
const pc = 3.08e18
const kpc = (10.0^3) * pc
const G = 6.67e-8
const gamma = 5.0 / 3.0

rvir=0.88*1.e3*kpc
r500 = 5*rvir/7
r_max=r500

include(string(main, "read.jl"))
include(string(main, "vel_module3D.jl"))
include(string(main, "nfluct.jl"))

d, d_dm, temp, vx, vy, vz, Npoint, conv_d  = read(snap, cluster)
#v = vmodule3D(vx,vy,vz,Npoint)

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

Ncell = floor(Int32, (r500/(20*kpc)))
include(string(main, "z_Projection.jl"))
include(string(main, "vel_module2D.jl"))
dx=0.02

d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D = projection(d, temp, vx, vy, vz, Npoint)

#********************* Let's consider only ************************
#********************* the cells inside R500 **********************

#for a box filtering we need to consider more cells used to calculate the mean

Ncell = floor(Int32, (r500/(20*kpc)))
Ncell_box= 15
Nhalf_box=floor(Int32,Ncell_box/2)

d = d[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
      yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box,
      zc-Ncell-Nhalf_box:zc+Ncell+Nhalf_box]

temp = temp[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
            yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box,
            zc-Ncell-Nhalf_box:zc+Ncell+Nhalf_box]

d_dm = d_dm[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
            yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box,
            zc-Ncell-Nhalf_box:zc+Ncell+Nhalf_box]


d_2D = d_2D[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
            yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box]

temp_2D = spec_T_2D[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
                    yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box]

vx_2D = vx_2D[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
              yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box]

vy_2D = vy_2D[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
           yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box]

vz_2D = vz_2D[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
           yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box]

xc = floor(Int32, length(d_2D[:,1])/2 +1)
yc = floor(Int32, length(d_2D[1,:])/2 +1)

include(string(main, "radial_profile.jl"))

Npoint=length(d_2D[1,:])
Nhalf = Int(floor(Npoint/2))
r_max = 20. * kpc * Ncell
#r_max=r500
rmin = 20. * kpc
radial_bin = Int(floor((1/0.1)*r_max/r500 + 1))

include(string(main, "2D_radial_profile.jl"))
d_2D_prof, temp_2D_prof, vx_2D_prof, vy_2D_prof, vz_2D_prof, radius = profile2D(d_2D, temp_2D, vx_2D, vy_2D, vz_2D, xc, yc, zc, r500, Nhalf, radial_bin, 0.1)

#****************** FILTERING *********************************

include(string(main, "2D_radial_filtering.jl"))
include(string(main, "v_mean2D.jl"))
deltad_2D, deltaT_2D, deltav_2D = filtering2D(d_2D, d_2D_prof, spec_T_2D, temp_2D_prof,
                                              vx_2D, vx_2D_prof, vy_2D, vy_2D_prof, vz_2D, vz_2D_prof,
                                              radial_bin, Npoint, xc,yc,zc, r500, 0.1)


v_2D = zeros(length(vx_2D[:,1]), length(vx_2D[:,1]))
for i in 1:length(vx_2D[:,1])
    for j in 1:length(vx_2D[:,1])

          v_2D[i,j] = sqrt(vx_2D[i,j]^2 + vy_2D[i,j]^2 + vz_2D[i,j]^2)

    end
end
#d_max, index = findmax(d_dm[30:40,30:40])
#xc = 30+index[1]
#yc = 30+index[2]

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

Npoint=length(d[1,1,:])

dx=0.02
subplot(1,2,1)
pcolormesh((-(Ncell+Nhalf_box)*dx:dx:(Ncell+Nhalf_box)*dx),(-(Ncell+Nhalf_box)*dx:dx:(Ncell+Nhalf_box)*dx),(((v_2D[:,:]))))

axis("scaled")
#title("IT90_3", fontsize= "20")
xticks(fontsize= 20)
yticks(fontsize=20)
cbar=colorbar(orientation="vertical")
cbar.ax.tick_params(labelsize=20)
cbar.ax.yaxis.get_offset_text().set(size=20)
#cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"Log \rho",fontsize= 25)
xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 30)
ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 30)

subplot(1,2,2)
pcolormesh((-(Ncell+Nhalf_box)*dx:dx:(Ncell+Nhalf_box)*dx),(-(Ncell+Nhalf_box)*dx:dx:(Ncell+Nhalf_box)*dx),(((deltav_2D[:,:]))))

axis("scaled")
#title("IT90_3", fontsize= "20")
xticks(fontsize= 20)
yticks(fontsize=20)
cbar=colorbar(orientation="vertical")
cbar.ax.tick_params(labelsize=20)
cbar.ax.yaxis.get_offset_text().set(size=20)
#cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"Log T",fontsize= 25)
xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 30)
ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 30)
