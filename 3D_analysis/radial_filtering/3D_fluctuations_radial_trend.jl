using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit

main = "/home/marco/Scrivania/Ettori_project/Codes/3D_analysis/radial_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"

name = ["IT90_0", "IT90_1", "IT90_2", "IT90_3", "IT90_4", "IT92_0", "IT92_1", "IT92_2"]
snap = ["0187", "0196", "0195", "0193", "0199", "0243", "0228", "0242"]
r_vir = [0.88, 1.29, 0.99, 0.86, 0.78, 1.42, 0.96, 1.01 ]

const mu = 0.62
const mp = 1.67e-24
const kb = 1.28e-16
const pc = 3.08e18
const kpc = (10.0^3) * pc
const G = 6.67e-8
const gamma = 5.0 / 3.0

include(string(main, "read.jl"))
include(string(main, "vel_module3D.jl"))
include(string(main, "nfluct.jl"))

for l in 1:length(snap)

rvir=r_vir[l]*1.e3*kpc
r500 = 5*rvir/7
thickness = 0.2
dr = thickness*r500 #thick of shell for radial profile CCHANGE THE FOLDER TOO IF YOU CHANGE THE THICKNESS
#r_max = 160 * 20 *kpc * sqrt(3)
r_max=r500


d, d_dm, temp, vx, vy, vz, Npoint, conv_d  = read(snap[l], cluster)
#v = vmodule3D(vx,vy,vz,Npoint)

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

N = thickness*200
Ncell = floor(Int32, (r500/(20*kpc))+ N)
#********************* Let's consider only ************************
#********************* the cells inside R500 **********************
d_dm = d_dm[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
d = d[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
temp = temp[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
vx = vx[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
vy = vy[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
vz = vz[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

Npoint=length(d[1,1,:])

dx=0.02
#subplot(1,2,1)
#pcolormesh((-(Ncell)*dx:dx:(Ncell)*dx),(-(Ncell)*dx:dx:(Ncell)*dx),(((d[:,:,zc]))))

#axis("scaled")
#title("IT90_3", fontsize= "20")
#xticks(fontsize= 20)
#yticks(fontsize=20)
#cbar=colorbar(orientation="vertical")
#cbar.ax.tick_params(labelsize=20)
#cbar.ax.yaxis.get_offset_text().set(size=20)
#cbar[:set_yticklabel](labelsize="large")
#cbar[:set_label](L"Log T",fontsize= 25)
#xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 30)
#ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 30)

#****************** RADIAL PROFILE **println("r=",r)**************************
include(string(main, "radial_profile.jl"))

Nhalf = Int(floor(Npoint/2))
r_max = 20. * kpc * Ncell * sqrt(3)
#r_max=r500
rmin = 20. * kpc
radial_bin = Int(floor(10*r_max/r500))
d_prof, temp_prof,vx_prof, vy_prof, vz_prof, radius =profile(d, temp, vx, vy, vz, xc, yc, zc, r500, rmin, Nhalf,radial_bin, thickness)##

bin=[1,2,3,4,5,6,7,8,9,10]
#output1=string("/home/marco/Scrivania/Ettori_project/", name[l],"/01_radial_filtering/", name[l],"_3D_temperature_radial_profile.txt")
#out1=open(output1,"w") do out1
#writedlm( out1, [bin temp_prof[1:10]])
#end

#****************** FILTERING *********************************

include(string(main, "radial_filtering.jl"))
deltav, deltad, deltaT= v_turb2(d, d_prof, vx, vx_prof, vy, vy_prof, vz, vz_prof, temp, temp_prof, Npoint, radial_bin, xc,yc,zc, thickness, r500)
println("3D filtering done")

rmax=r500
Ncell = floor(Int32, (r500/(20*kpc)))
radial_bin = Int(floor(10*rmax/r500))
d_dm = d_dm[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
d = d[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
temp = temp[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
vx = vx[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
vy = vy[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
vz = vz[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
deltad = deltad[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
deltaT = deltaT[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
deltav = deltav[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]

xc = floor(Int32, length(d_dm[:,1,1])/2 +1)
yc = floor(Int32, length(d_dm[1,1,:])/2 +1)
zc = floor(Int32, length(d_dm[1,1,:])/2 +1)

#subplot(1,2,2)
#pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((deltad[:,:,zc-N]))))

#axis("scaled")
#title("IT90_3", fontsize= "20")
#xticks(fontsize= 20)
#yticks(fontsize=20);
#cbar=colorbar(orientation="vertical")
#cbar.ax.tick_params(labelsize=20)
#cbar.ax.yaxis.get_offset_text().set(size=20)
#cbar[:set_yticklabel](labelsize="large")
#cbar[:set_label](L"Log(\delta T)",fontsize= 25)
#xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 30)
#ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 30)

#readline()
#clf()

rmax=r500
radial_bin = Int(floor(10*rmax/r500))

#include(string(main, "3D_average_radial_perturbation.jl"))
#dens_shell_mean, temp_shell_mean, v_shell_mean, v_shell_mean2=average_radial_perturb(d, d_prof,temp, temp_prof,
#                                              vx, vx_prof, vy, vy_prof, vz, vz_prof, radial_bin, xc, yc ,zc, r500, thickness)

include(string(main, "3D_slope_radial_profile.jl"))
slope_dt, slope_dv, slope_tv, slope_tv2=slope_radial_profile(deltad,deltaT,deltav, radial_bin, xc, yc ,zc, r500, thickness)

bin=[1,2,3,4,5,6,7,8,9,10]
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/01_radial_filtering/new_",name[l],"_3D_radial_filtering_fluctuations_radial_trend.txt")
#out1=open(output1,"w") do out1
#writedlm( out1, [bin dens_shell_mean temp_shell_mean v_shell_mean v_shell_mean2])
#end

output1=file1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/01_radial_filtering/new_",name[l],"_3D_radial_filtering_slope_radial_trend.txt")
out1=open(output1,"w") do out1
writedlm( out1, [bin slope_dt slope_dv slope_tv slope_tv2])
end

end
