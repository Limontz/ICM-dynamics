using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit
using HypothesisTests

main = "/home/marco/Scrivania/Ettori_project/Codes/3D_analysis/radial_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"

include(string(main, "read.jl"))
include(string(main, "vel_module3D.jl"))
include(string(main, "nfluct.jl"))
include(string(main, "radial_profile.jl"))
include(string(main, "radial_filtering.jl"))


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

sigmaT = zeros(length(snap))
sigman = zeros(length(snap))

for l in 1:length(snap)

rvir=r_vir[l]*1.e3*kpc
r500 = 5*rvir/7
thickness = 0.2
dr = thickness*r500 #thick of shell for radial profile
#r_max = 160 * 20 *kpc * sqrt(3)
r_max=r500

d, d_dm, temp, vx, vy, vz, Npoint, conv_d  = read(snap[l], cluster)
#v = vmodule3D(vx,vy,vz,Npoint)

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

N = floor(Int32, thickness*200)
Ncell = floor(Int32, (r500/(20*kpc)) + N)
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


dmin, index2 = findmin(d[:,:,zc])
dmax, index3 = findmax(d[:,:,zc])

Npoint=length(d[1,1,:])

#************************* CLUMPS FILTERING *********************************
#include(string(main, "clump_excision.jl"))
#include(string(main, "clump_masking.jl"))
#include(string(main, "pdf.jl"))


#radial_bin = Int(floor(10*r_max/r500))
#lim95, lim99, shell_median = masking(d, Npoint, xc, yc, zc, radial_bin)
#d = excision(d, temp, lim99, shell_median, xc, yc, zc, Npoint)

#****************** RADIAL PROFILE **println("r=",r)**************************


Nhalf = Int(floor(Npoint/2))
r_max = 20. * kpc * Ncell
#r_max=r500
rmin = 20. * kpc
radial_bin = Int(floor((1/thickness)*r_max/r500 + 1))

d_prof, temp_prof,vx_prof, vy_prof, vz_prof, radius =profile(d, temp, vx, vy, vz, xc, yc, zc, r500, rmin, Nhalf,radial_bin, thickness)##

#****************** FILTERING *********************************

dx=0.02
#subplot(1,2,1)
#pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),((log10.(d[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc]))), vmin=log10(dmin), vmax=log10(dmax))

#axis("scaled")
#title("IT90_3", fontsize= "20")
#xticks(fontsize= 20)
#yticks(fontsize=20)
#cbar=colorbar(orientation="vertical")
#cbar.ax.tick_params(labelsize=20)
#cbar.ax.yaxis.get_offset_text().set(size=20)
##cbar[:set_yticklabel](labelsize="large")
#cbar[:set_label](L"Log T",fontsize= 25)
#xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 30)
#ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 30)

deltav, deltad, deltaT= v_turb2(d, d_prof, vx, vx_prof, vy, vy_prof, vz, vz_prof, temp, temp_prof, Npoint, radial_bin, xc,yc,zc, thickness,r500)
println("3D filtering done")

Ncell = floor(Int32, (r500/(20*kpc)))


d_dm = d_dm[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
deltad = deltad[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
deltaT = deltaT[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]
deltav = deltav[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell,zc-Ncell:zc+Ncell]

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

#sigmaT[l] = std(deltaT)
#sigman[l] = std(deltad)

#end

#output1=file1=string("/home/marco/Scrivania/Ettori_project/Plots/3D_02_radial_filtering_mean_fluctuations.txt")
#out1=open(output1,"w") do out1
#writedlm( out1, [sigmaT sigman])
#end

#subplot(1,2,2)
#pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((deltad[:,:,zc]))))

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

#************************ masking the largest fluctuations ***********************

#nr_bin = 100
#n, bin, max = pdf(deltad, nr_bin)
#plot(bin, n)
#r_max=r500
#radial_bin = Int(floor(10*r_max/r500))
#lim95, lim99, shell_median = masking(deltad, Npoint, Ncell+1, Ncell+1, Ncell+1, radial_bin, r_max)
#print(lim99)
#deltad, deltaT, deltav = mask(deltad, deltaT, deltav, Ncell+1, Ncell+1, Ncell+1, lim99, radial_bin, r_max)

#nr_bin = 100
#n, bin, max = pdf(deltad, nr_bin)
#plot(bin, n)
#readline()

#************************* LINEAR FIT ********************************************


deltaT = vec(deltaT)
deltav = vec(deltav)
deltad = vec(deltad)

println("vectorization done")

println("----------- 3D FLUCTUATIONS: temperature-density -------------")

m(x, p) = p[1] .* x

p0=[1.]
fit = curve_fit(m, ((deltad)), ((deltaT)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid

x = collect(minimum((deltad)):0.001:maximum((deltad)))
y = zeros(length(x))

for i in 1:length(x)
    y[i] = param[1]*((x[i]))
    #y[i] = param[1]*(abs(x[i]))
end

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
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/02_radial_filtering/new_",name[l],"_3D_radial_filtering_fit_parameters.txt")
#out1=open(output1,"w") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1=string("/home/marco/Scrivania/Ettori_project/02_3D_radial_filtering/3D_02_radial_filtering_density_temperature_fit_parameters.txt")
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end

#readline()
#clf()

println("----------- 3D FLUCTUATIONS: density-velocity-------------")

m(x, p) = p[1].* x
p0=[1.]
fit = curve_fit(m, (abs.(deltad)), (abs.(deltav)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid


x = collect(minimum(abs.(deltad)):0.001:maximum(abs.(deltad)))
y = zeros(length(x))

for i in 1:length(x)
    #y[i] = param[1] + param[2]*(abs(x[i]))
    y[i] = param[1]*(abs(x[i]))
end

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

R = Statistics.cor(abs.(deltad), abs.(deltav))
p_value=pvalue(CorrelationTest(abs.(deltad),abs.(deltav)))

#ans = "3D FLUCTUATIONS density-velocity:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/02_radial_filtering/new_",name[l],"_3D_radial_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1=string("/home/marco/Scrivania/Ettori_project/02_3D_radial_filtering/3D_02_radial_filtering_density_velocity_fit_parameters.txt")
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end



#readline()
#clf()

println("----------- 3D FLUCTUATIONS: temperature-velocity -------------")

m(x, p) = p[1].* x
p0=[1.]
fit = curve_fit(m, (abs.(deltaT)), (abs.(deltav)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid


x = collect(minimum(abs.(deltaT)):0.001:maximum(abs.(deltaT)))
y = zeros(length(x))

for i in 1:length(x)
    #y[i] = param[1] + param[2]*(abs(x[i]))
    y[i] = param[1]*(abs(x[i]))
end

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
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/02_radial_filtering/new_",name[l],"_3D_radial_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1=string("/home/marco/Scrivania/Ettori_project/02_3D_radial_filtering/3D_02_radial_filtering_temperature_velocity_fit_parameters.txt")
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end

#readline()
#clf()

println("----------- 3D FLUCTUATIONS: temperature-velocity^2 -------------")

m(x, p) = p[1] .* x
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

R = Statistics.cor(abs.(deltaT), abs.(deltav).^2)
p_value=pvalue(CorrelationTest(abs.(deltaT),abs.(deltav).^2))

#ans = "3D FLUCTUATIONS temperature-velocity^2:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/02_radial_filtering/new_",name[l],"_3D_radial_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1="/home/marco/Scrivania/Ettori_project/02_3D_radial_filtering/3D_02_radial_temperature_velocity2_filtering_fit_parameters.txt"
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end

end
