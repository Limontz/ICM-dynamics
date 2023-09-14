using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit


main = "/home/marco/Scrivania/Ettori_project/Codes/3D_analysis/box_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"



const mu = 0.62
const mp = 1.67e-24
const kb = 1.28e-16
const pc = 3.08e18
const kpc = (10.0^3) * pc
const G = 6.67e-8
const gamma = 5.0 / 3.0

name = ["IT90_0", "IT90_1", "IT90_2", "IT90_3", "IT90_4", "IT92_0", "IT92_1", "IT92_2"]
snap = ["0187", "0196", "0195", "0193", "0199", "0243", "0228", "0242"]
r_vir = [0.88, 1.29, 0.99, 0.86, 0.78, 1.42, 0.96, 1.01 ]

include(string(main, "read.jl"))
include(string(main, "vel_module3D.jl"))
include(string(main, "nfluct.jl"))

for l in 1:length(snap)

rvir=r_vir[l]*1.e3*kpc
r500 = 5*rvir/7
#dr = 0.1*r500 #thick of shell for radial profile
#r_max = 160 * 20 *kpc * sqrt(3)
r_max=r500


d, d_dm, temp, vx, vy, vz, Npoint, conv_d  = read(snap[l], cluster)
#v = vmodule3D(vx,vy,vz,Npoint)

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

Ncell = floor(Int32, (r500/(20*kpc)))

#********************* Let's consider only ************************
#********************* the cells inside R500 **********************

#for a box filtering we need to consider more cells used to calculate the mean
Ncell_box= 15
Nhalf_box=floor(Int32,Ncell_box/2)

println("here1")
d_dm = d_dm[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
      yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box,
      zc-Ncell-Nhalf_box:zc+Ncell+Nhalf_box]

d = d[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
      yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box,
      zc-Ncell-Nhalf_box:zc+Ncell+Nhalf_box]

temp = temp[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
            yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box,
            zc-Ncell-Nhalf_box:zc+Ncell+Nhalf_box]

vx = vx[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
        yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box,
        zc-Ncell-Nhalf_box:zc+Ncell+Nhalf_box]

vy = vy[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
        yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box,
        zc-Ncell-Nhalf_box:zc+Ncell+Nhalf_box]

vz = vz[xc-Ncell-Nhalf_box:xc+Ncell+Nhalf_box,
        yc-Ncell-Nhalf_box:yc+Ncell+Nhalf_box,
        zc-Ncell-Nhalf_box:zc+Ncell+Nhalf_box]
println("here2")
d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

Npoint=length(d[1,1,:])

#dx=0.02
#subplot(1,2,1)
#pcolormesh((-(Ncell+Nhalf_box)*dx:dx:(Ncell+Nhalf_box)*dx),(-(Ncell+Nhalf_box)*dx:dx:(Ncell+Nhalf_box)*dx),(((d[:,:,zc]))))

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


#****************** FILTERING *********************************
include(string(main, "v_turb2.jl"))
include(string(main, "v_mean2.jl"))
deltav, deltad, deltaT, L = v_turb2(d, vx, vy, vz, temp, Npoint, Nhalf_box, xc,yc,zc)
println("3D filtering done")

Ncell = floor(Int32, (r500/(20*kpc)))

println("here3")
d_dm = d_dm[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]

deltaT = deltaT[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
deltav = deltav[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
deltad = deltad[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
println("here4")
d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

#subplot(1,2,2)
#pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((deltad[:,:,zc-Nhalf_box]))))

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

include(string(main, "3D_average_radial_perturbation.jl"))
dens_shell_mean, temp_shell_mean, v_shell_mean, slope_dt, slope_dv, slope_tv, slope_tv2=average_radial_perturb(deltad, deltaT,deltav, radial_bin, xc, yc ,zc, r500)


bin=[1,2,3,4,5,6,7,8,9,10]
#output1=string("/home/marco/Scrivania/Ettori_project/IT92_2/300_kpc_filtering/IT92_2_3D_box_filtering_fluctuations_radial_trend.txt")
#out1=open(output1,"a+") do out1
#writedlm( out1, [bin dens_shell_mean temp_shell_mean v_shell_mean])
#end

output1=file1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/300_kpc_filtering/new_",name[l],"_3D_box_filtering_slope_radial_trend.txt")
out1=open(output1,"w") do out1
writedlm( out1, [bin slope_dt slope_dv slope_tv slope_tv2])
end

end
