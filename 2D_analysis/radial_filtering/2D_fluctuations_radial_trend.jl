using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit


main = "/home/marco/Scrivania/Ettori_project/Codes/2D_analysis/radial_filtering/"
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

for l in 1:length(name)

rvir=r_vir[l]*1.e3*kpc
r500 = 5*rvir/7
thickness = 0.2
dr = thickness*r500 #thick of shell for radial profile  CHANGE THE FOLDER TOO IF YOU CHANGE THE THICKNESS
#r_max = 160 * 20 *kpc * sqrt(3)
r_max=r500

d, d_dm, temp, vx, vy, vz, Npoint, conv_d  = read(snap[l], cluster)
#v = vmodule3D(vx,vy,vz,Npoint)

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]


#********************* Let's consider only ************************
#********************* the cells inside R500 **********************
N = thickness*200
Ncell = floor(Int32, (r500/(20*kpc)) + N)

d_dm = d_dm[xc-Ncell:xc+Ncell,
            yc-Ncell:yc+Ncell,
            zc-Ncell:zc+Ncell]


d = d[xc-Ncell:xc+Ncell,
      yc-Ncell:yc+Ncell,
      zc-Ncell:zc+Ncell]

temp = temp[xc-Ncell:xc+Ncell,
            yc-Ncell:yc+Ncell,
            zc-Ncell:zc+Ncell]

vx = vx[xc-Ncell:xc+Ncell,
            yc-Ncell:yc+Ncell,
            zc-Ncell:zc+Ncell]

vy = vy[xc-Ncell:xc+Ncell,
            yc-Ncell:yc+Ncell,
            zc-Ncell:zc+Ncell]

vz = vz[xc-Ncell:xc+Ncell,
            yc-Ncell:yc+Ncell,
            zc-Ncell:zc+Ncell]

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

Npoint=length(d[:,1,1])


#************************* CLUMPS FILTERING *********************************
#include(string(main, "clump_excision.jl"))
#include(string(main, "clump_masking.jl"))
#include(string(main, "pdf.jl"))


#radial_bin = Int(floor(10*r_max/r500))
#lim95, lim99, shell_median = masking(d, Npoint, xc, yc, zc, radial_bin)
#d = excision(d, temp, lim99, shell_median, xc, yc, zc, Npoint)


#*********************PROJECTION***********************************

include(string(main, "y_Projection.jl"))
include(string(main, "vel_module2D.jl"))
dx=0.02
d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D = projection(d, temp, vx, vy, vz, Npoint)

dx=0.02
#subplot(1,2,1)
#pcolormesh(((-Ncell)*dx:dx:(Ncell)*dx),((-Ncell)*dx:dx:(Ncell)*dx),(((d_2D))))

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

#****************** RADIAL PROFILE ****************************

Npoint=length(d_2D[1,:])
Nhalf = Ncell
#r_max = 20. * kpc * Ncell * sqrt(3)
#r_max=r500
rmin = 20. * kpc
radial_bin = 10



include(string(main, "2D_radial_profile.jl"))
d_2D_prof, temp_2D_prof, vx_2D_prof, vy_2D_prof, vz_2D_prof, radius = profile2D(d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D, xc, yc, zc, r500, Nhalf, radial_bin, thickness)

#****************** FILTERING *********************************thickness

include(string(main, "2D_radial_filtering.jl"))
include(string(main, "v_mean2D.jl"))
deltad_2D, deltaT_2D, deltav_2D = filtering2D(d_2D, d_2D_prof, spec_T_2D, temp_2D_prof,
                                              vx_2D, vx_2D_prof, vy_2D, vy_2D_prof, vz_2D, vz_2D_prof,
                                              radial_bin, Npoint, xc,yc,zc, r500, thickness)


println("2D filtering done")

Ncell = floor(Int32, (r500/(20*kpc)))

deltad_2D = deltad_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
deltaT_2D = deltaT_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
deltav_2D = deltav_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]

d_2D = d_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
spec_T_2D = spec_T_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
vx_2D = vx_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
vy_2D = vy_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
vz_2D = vz_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]


xc = floor(Int32, length(d_2D[:,1])/2 +1)
yc = floor(Int32, length(d_2D[1,:])/2 +1)



#subplot(1,2,2)
#pcolormesh(((-Ncell)*dx:dx:(Ncell)*dx),((-Ncell)*dx:dx:(Ncell)*dx),(((deltad_2D))))

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


#include(string(main, "2D_average_radial_perturbation.jl"))
#dens_shell_mean, dens_shell_mean2, temp_shell_mean, v_shell_mean, v_shell_mean2=average_radial_perturb(d_2D, d_2D_prof, spec_T_2D, temp_2D_prof,
                                              #vx_2D, vx_2D_prof, vy_2D, vy_2D_prof, vz_2D, vz_2D_prof, radial_bin, xc, yc ,zc, r500, thickness)

include(string(main, "2D_slope_radial_profile.jl"))
slope_dt, slope_dv, slope_tv, slope_tv2=slope_radial_profile(deltad_2D,deltaT_2D,deltav_2D, radial_bin, xc, yc ,zc, r500, thickness)


bin=[1,2,3,4,5,6,7,8,9,10]
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/01_radial_filtering/new_y_",name[l],"_2D_radial_filtering_fluctuations_radial_trend.txt")
#out1=open(output1,"w") do out1
#writedlm( out1, [bin dens_shell_mean dens_shell_mean2 temp_shell_mean v_shell_mean v_shell_mean2])
#end



output1=string("/home/marco/Scrivania/Ettori_project/", name[l],"/02_radial_filtering/y_", name[l],"_2D_radial_filtering_slope_radial_trend.txt")
out1=open(output1,"w") do out1
writedlm( out1, [bin slope_dt slope_dv slope_tv slope_tv2])
end

end
