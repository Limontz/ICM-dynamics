using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit
using HypothesisTests

name = ["IT90_0", "IT90_1", "IT90_2", "IT90_3", "IT90_4", "IT92_0", "IT92_1", "IT92_2"]
snap = ["0187", "0196", "0195", "0193", "0199", "0243", "0228", "0242"]
r_vir = [0.88, 1.29, 0.99, 0.86, 0.78, 1.42, 0.96, 1.01 ]

main = "/home/marco/Scrivania/Ettori_project/Codes/3D_analysis/box_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"

include(string(main, "read.jl"))
include(string(main, "vel_module3D.jl"))
include(string(main, "nfluct.jl"))
include(string(main, "v_turb2.jl"))
include(string(main, "v_mean2.jl"))

const mu = 0.62
const mp = 1.67e-24
const kb = 1.28e-16
const pc = 3.08e18
const kpc = (10.0^3) * pc
const G = 6.67e-8
const gamma = 5.0 / 3.0

sigmaT = zeros(length(snap))
sigman = zeros(length(snap))

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

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

Npoint=length(d[1,1,:])

dx=0.02
#subplot(1,2,1)
#pcolormesh((-(Ncell+Nhalf_box)*dx:dx:(Ncell+Nhalf_box)*dx),(-(Ncell+Nhalf_box)*dx:dx:(Ncell+Nhalf_box)*dx),(((sqrt(3) .* vx[:,:,zc]))))

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

#************************* CLUMPS FILTERING *********************************
#include(string(main, "clump_excision.jl"))
#include(string(main, "clump_masking.jl"))
#include(string(main, "pdf.jl"))


#radial_bin = Int(floor(10*r_max/r500))
#lim95, lim99, shprintln(Nhalf_box)ell_median = masking(d, Npoint, xc, yc, zc, radial_bin)
#d = excision(d, temp, lim99, shell_median, xc, yc, zc, Npoint)


#****************** RADIAL PROFILE ****************************
#include(string(main, "radial_profile.jl"))

#Nhalf = Int(floor(Npoint/2))
#r_max = 20. * kpc * Nhalf * sqrt(3)
#r_max=r500
#rmin = 20. * kpc
#radial_bin = Int(floor(10*r_max/r500))
#d_prof, temp_prof, radius =profile(d, temp, xc, yc, zc, r500, rmin, Nhalf,radial_bin)##

#****************** FILTERING *********************************
print("here")
deltav, deltad, deltaT, L = v_turb2(d, vx, vy, vz, temp, Npoint, Nhalf_box, xc,yc,zc)
println("3D filtering done")

Ncell = floor(Int32, (r500/(20*kpc)))
d_dm = (d_dm[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell])
deltaT = (deltaT[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell])
deltav = (deltav[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell])
deltad = (deltad[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell])

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

sigmaT[l] = std(deltaT)
sigman[l] = std(deltad)

end

output1=file1=string("/home/marco/Scrivania/Ettori_project/Plots/3D_300_box_filtering_mean_fluctuations.txt")
out1=open(output1,"w") do out1
writedlm( out1, [sigmaT sigman])
end


dx=0.02

#subplot(1,2,2)
#pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((deltaT[:,:,zc-Nhalf_box]))))

#axis("scaled")
##title("IT90_3", fontsize= "20")
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

#************************ masking the largest fluctuations ***********************

#include(string(main, "clump_excision.jl"))
#include(string(main, "clump_masking.jl"))
#include(string(main, "pdf.jl"))
#include(string(main, "large_fluctuations_mask_2D.jl"))


#radial_bin = Int(floor(10*r_max/r500))
#lim95, lim99, shell_median = masking(deltad, Npoint, xc, yc, zc, radial_bin)
#deltad, deltaT = mask2D(deltad, deltaT, xc,yc, zc, lim99, Npoint, radial_bin)

#-------------- Fluctuations as result of filtering --------------------------------

deltaT = vec(deltaT)
deltav = vec(deltav)
deltad = vec(deltad)

println("vectorization done")
#************************* LINEAR FIT ********************************************

println("----------- 3D FLUCTUATIONS: temperature-density -------------")

m(x, p) = p[1] .* x

p0=[1.]
fit = curve_fit(m, ((deltad)), ((deltaT)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid


#x = collect(minimum((deltad)):0.001:maximum((deltad)))
#y = zeros(length(x))

#for i in 1:length(x)
 #   y[i] = param[1]*((x[i]))
    #y[i] = param[1]*(abs(x[i]))
#end

#hist, xedges, yedges=hist2D(((deltad)),((deltaT)), bins=(10 .^(-3:0.2:0), 10 .^(-3:.2:0)), cmap="afmhot_r")
#hist, xedges, yedges=hist2D(((deltad)),((deltaT)), bins=(40,40), cmap="afmhot_r")
#plot(deltad_2D, deltaT_2D, linestyle="", marker=".")
#plot((x), (y))
#cbar=colorbar(orientation="vertical")
#ylabel(L"(\delta T / <T>)_{3D}", fontsize=15)
#xlabel(L"(\delta \rho / <\rho>)_{3D}", fontsize=15)
#xticks(fontsize=15)
#yticks(fontsize=15)
#yscale("log")
#xscale("log")
#xlim(-0.3,0.3)
#ylim(-0.3,0.3)

R = Statistics.cor(deltad, deltaT)
p_value=pvalue(CorrelationTest(deltad,deltaT))

#ans = "3D FLUCTUATIONS density-temperature:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/100_kpc_filtering/new_",name[l],"_3D_box_filtering_fit_parameters.txt")
#out1=open(output1,"w") do out1
#write( out1, ans)
#end
#close(output1)

#output1=file1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/100_kpc_filtering/new_",name[l],"_3D_box_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#writedlm( out1, [param[1] sigma_par[1] R p_value])
#end


println("----------- 3D FLUCTUATIONS: density-velocity-------------")

m(x, p) = p[1].* x
p0=[1.]
fit = curve_fit(m, (abs.(deltad)), (abs.(deltav)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid

R = Statistics.cor(abs.(deltad), abs.(deltav))
p_value=pvalue(CorrelationTest(abs.(deltad),abs.(deltav)))

#ans = "3D FLUCTUATIONS density-velocity:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/100_kpc_filtering/new_",name[l],"_3D_box_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

#output1=file1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/100_kpc_filtering/new_",name[l],"_3D_box_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#writedlm( out1, [param[1] sigma_par[1] R p_value])
#end


#x = collect(minimum(abs.(deltad)):0.001:maximum(abs.(deltad)))
#y = zeros(length(x))

#for i in 1:length(x)
 #   y[i] = param[1]*(abs(x[i]))
    #y[i] = param[1]*(abs(x[i]))
#end

#hist, xedges, yedges=hist2D((abs.(deltad)),(abs.(deltav)), bins=(10 .^(-3:0.2:0.5), 10 .^(-3:.2:0.5)), cmap="afmhot_r")
#plot((x), (y))
#cbar=colorbar(orientation="vertical")
#ylabel(L"(\delta v / <v>)_{3D}", fontsize=15)
#xlabel(L"(\delta \rho / <\rho>)_{3D}", fontsize=15)
#xticks(fontsize=15)
#yticks(fontsize=15)
#yscale("log")
#xscale("log")
#xlim(-3,0)
#ylim(-3,0)

#clf()

println("----------- 3D FLUCTUATIONS: temperature-velocity -------------")

p0=[1.]
fit = curve_fit(m, (abs.(deltaT)), (abs.(deltav)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid

#x = collect(minimum(abs.(deltaT)):0.001:maximum(abs.(deltaT)))
#y = zeros(length(x))

#for i in 1:length(x)
 #   y[i] = param[1]*(abs(x[i]))
    #y[i] = param[1]*(abs(x[i]))
#end

#hist, xedges, yedges=hist2D((abs.(deltaT)),(abs.(deltav)), bins=(10 .^(-3:0.2:0.5), 10 .^(-3:.2:0.5)), cmap="afmhot_r")
#plot((x), (y))
#cbar=colorbar(orientation="vertical")
#ylabel(L"(\delta v / <v>)_{3D}", fontsize=15)
#xlabel(L"(\delta T / <T>)_{3D}", fontsize=15)
#xticks(fontsize=15)
#yticks(fontsize=15)
#yscale("log")
#xscale("log")
#xlim(-3,0)
#ylim(-3,0)

R = Statistics.cor(abs.(deltaT), abs.(deltav))
p_value=pvalue(CorrelationTest(abs.(deltaT),abs.(deltav)))

#ans = "3D FLUCTUATIONS temperature-velocity:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/100_kpc_filtering/new_",name[l],"_3D_box_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

#output1=file1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/100_kpc_filtering/new_",name[l],"_3D_box_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#writedlm( out1, [param[1] sigma_par[1] R p_value])
#end

println("----------- 3D FLUCTUATIONS: temperature-velocity^2 -------------")

p0=[1.]
fit = curve_fit(m, (abs.(deltaT)), (abs.(deltav)).^2, p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid

println(param)

x = collect(minimum(abs.(deltaT)):0.001:maximum(abs.(deltaT)))
y = zeros(length(x))

for i in 1:length(x)
   y[i] = param[1]*(abs(x[i]))
    #y[i] = param[1]*(abs(x[i]))
end

#hist, xedges, yedges=hist2D((abs.(deltaT)),(abs.(deltav)).^2, bins=(10 .^(-3:0.2:0.5), 10 .^(-3:.2:0.5)), cmap="afmhot_r")
#plot((x), (y), linestyle="-", color="blue")
#cbar=colorbar(orientation="vertical")
#ylabel(L"(\delta v / <v>)_{3D}", fontsize=15)
#xlabel(L"(\delta T / <T>)_{3D}", fontsize=15)
#xticks(fontsize=15)
#yticks(fontsize=15)
#yscale("log")
#xscale("log")
#xlim(-3,0)
#ylim(-3,0)

R = Statistics.cor(abs.(deltaT), abs.(deltav).^2)
p_value=pvalue(CorrelationTest(abs.(deltaT),abs.(deltav).^2))


#ans = "3D FLUCTUATIONS temperature-velocity^2:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/100_kpc_filtering/new_",name[l],"_3D_box_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

#output1=file1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/100_kpc_filtering/new_",name[l],"_3D_box_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#writedlm( out1, [param[1] sigma_par[1] R p_value])
#end

#end
