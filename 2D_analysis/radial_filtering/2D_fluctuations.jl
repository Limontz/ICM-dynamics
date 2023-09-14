using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit
using HypothesisTests


main = "/home/marco/Scrivania/Ettori_project/Codes/2D_analysis/radial_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"

include(string(main, "read.jl"))
include(string(main, "vel_module3D.jl"))
include(string(main, "nfluct.jl"))
include(string(main, "x_Projection.jl"))
include(string(main, "vel_module2D.jl"))
include(string(main, "radial_profile.jl"))
include(string(main, "2D_radial_profile.jl"))
include(string(main, "2D_radial_filtering.jl"))
include(string(main, "v_mean2D.jl"))

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
Ncell = floor(Int32, (r500/(20*kpc) + N))

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

#************************* CLUMPS FILTERING *********************************
#include(string(main, "clump_excision.jl"))
#include(string(main, "clump_masking.jl"))
#include(string(main, "pdf.jl"))


#radial_bin = Int(floor(10*r_max/r500))
#lim95, lim99, shell_median = masking(d, Npoint, xc, yc, zc, radial_bin)
#d = excision(d, temp, lim99, shell_median, xc, yc, zc, Npoint)


#*********************PROJECTION***********************************

dx=0.02
d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D = projection(d, temp, vx, vy, vz, Npoint)
#v_2D = vmodule2D(vx,vy,vz,Npoint)
#pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((d_2D))))

#axis("scaled")
#title("IT90_3", fontsize= "20")
#xticks(fontsize= 20)
#yticks(fontsize=20);
#cbar=colorbar(orientation="vertical")
#cbar.ax.tick_params(labelsize=20)
#cbar.ax.yaxis.get_offset_text().set(size=20)
#cbar[:set_yticklabel](labelsize="large")
#cbar[:set_label](L"Log(T [k])",fontsize= 25)
#xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 30)
#ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 30)
#readline()
#clf()
#****************** RADIAL PROFILE ****************************

Npoint=length(d_2D[1,:])
Nhalf = Int(floor(Npoint/2))
r_max = 20. * kpc * Ncell
#r_max=r500
rmin = 20. * kpc
radial_bin = Int(floor((1/thickness)*r_max/r500 + 1))

dx=0.02
#subplot(1,2,1)
#pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((d_2D))))
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


d_2D_prof, temp_2D_prof, vx_2D_prof, vy_2D_prof, vz_2D_prof, radius = profile2D(d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D, xc, yc, zc, r500, Nhalf, radial_bin, thickness)

#****************** FILTERING *********************************


deltad_2D, deltaT_2D, deltav_2D = filtering2D(d_2D, d_2D_prof, spec_T_2D, temp_2D_prof,
                                              vx_2D, vx_2D_prof, vy_2D, vy_2D_prof, vz_2D, vz_2D_prof,
                                              radial_bin, Npoint, xc,yc,zc, r500, thickness)


println("2D filtering done")

Ncell = floor(Int32, (r500/(20*kpc)))

deltad_2D = deltad_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
deltaT_2D = deltaT_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
deltav_2D = deltav_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]

xc = floor(Int32, length(d_2D[:,1])/2 +1)
yc = floor(Int32, length(d_2D[1,:])/2 +1)

sigmaT[l] = std(deltaT_2D)
sigman[l] = std(deltad_2D)


#output1=file1=string("/home/marco/Scrivania/Ettori_project/Plots/2D_02_radial_filtering_mean_fluctuations.txt")
#out1=open(output1,"w") do out1
#writedlm( out1, [sigmaT sigman])
#end

#exit()


#subplot(1,2,2)
#pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((deltad_2D))))

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
xc = floor(Int32, length(d_2D[:,1])/2 +1)
yc = floor(Int32, length(d_2D[1,:])/2 +1)

#readline()
#clf()

#************************ masking the largest fluctuations ***********************

#include(string(main, "clump_masking_2D.jl"))
#include(string(main, "pdf.jl"))
#include(string(main, "large_fluctuations_mask_2D.jl"))


radial_bin = Int(floor(10*r500/r500))
#lim95, lim99, shell_median = masking(deltad_2D, Npoint, xc, yc, zc, radial_bin, r500)

#deltad_2D, deltaT_2D, deltav_2D = mask2D(deltad_2D, deltaT_2D, deltav_2D, xc,yc, zc, lim95, Npoint, radial_bin, r500)

#nr_bin = 100
#n, bin, max = pdf(deltad_2D, nr_bin)
#plot(bin, n)
#readline()


#****************** DeltaT/T-deltav relation *********************************

Ncell = floor(Int32, (r500/(20*kpc)))
println((Ncell), " ", 0.6*r500/kpc)

#-------------- Fluctuations as result of filtering --------------------------------
#array_size=length(deltaT_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell])

deltaT_2D = vec(deltaT_2D)
deltad_2D = vec(deltad_2D)
deltav_2D = vec(deltav_2D)

println("vectorization done")
#************************* LINEAR FIT ********************************************

m(x, p) =  p[1] .* x
#m(x, p) = log.(abs.(1 .+ p[1] .* x[:,2].^4 .+ p[2] .* x[:,2].^2 .* x[:,1] .* x[:,2]) )
#m(x,p) = log.(1. .+ p[1]^2 .* x.^2)
p0=[1.]
fit = curve_fit(m, ((deltad_2D)), ((deltaT_2D)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid



x = collect(minimum((deltad_2D)):0.001:maximum((deltad_2D)))
y = zeros(length(x))

for i in 1:length(x)
    #y[i] = param[1] + param[2]*log10(abs(x[i]))
    y[i] = param[1]*((x[i]))
end

#h, xedges, yedges=hist2D(log10.(abs.(deltad_2D)),log10.(abs.(deltaT_2D)), bins=(30,30), cmap="afmhot_r")
#hist, xedges, yedges=hist2D(((deltad_2D)),((deltaT_2D)), bins=(30,30), cmap="afmhot_r")
#hist=hist2D((abs.(deltad_2D)),(abs.(deltaT_2D)), bins=(xedges,yedges), cmap="afmhot_r")
#plot(deltad_2D, deltaT_2D, linestyle="", marker=".")
#plot((x), (y))

#cbar=colorbar(orientation="vertical")
#ylabel(L"(\delta T /<T>)_{2D}", fontsize=15)
#xlabel(L"(\delta \rho / <\rho>)_{2D}", fontsize=15)
#xticks(fontsize=15)
#yticks(fontsize=15)
#yscale("log")
#xscale("log")
#xlim(-0.3,0.3)
#ylim(-0.3,0.3)

R = Statistics.cor(deltad_2D, deltaT_2D)
p_value=pvalue(CorrelationTest(deltad_2D,deltaT_2D))

#ans = "2D FLUCTUATIONS density-temperature:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/02_radial_filtering/new_",name[l],"_2D_radial_filtering_fit_parameters.txt")
#out1=open(output1,"w") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1=string("/home/marco/Scrivania/Ettori_project/02_2D_radial_filtering/x_2D_02_radial_filtering_density_temperature_fit_parameters.txt")
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end

#readline()
#clf()


println("-----------  2D FLUCTUATIONS: density-velocity -------------")

m(x, p) =  p[1] .* x


p0=[1.]
fit = curve_fit(m, (abs.(deltad_2D)), (abs.(deltav_2D)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid

x = collect(minimum(abs.(deltad_2D)):0.001:maximum(abs.(deltad_2D)))
y = zeros(length(x))

for i in 1:length(x)
    #y[i] = param[1] + param[2]*log10(abs(x[i]))
    y[i] = param[1]*(abs(x[i]))
end

R = Statistics.cor(abs.(deltad_2D), abs.(deltav_2D))
p_value=pvalue(CorrelationTest(abs.(deltad_2D),abs.(deltav_2D)))

#ans = "2D FLUCTUATIONS density-velocity:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/02_radial_filtering/new_",name[l],"_2D_radial_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1=string("/home/marco/Scrivania/Ettori_project/02_2D_radial_filtering/x_2D_02_radial_filtering_density_velocity_fit_parameters.txt")
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end



#hist, xedges, yedges=hist2D((abs.(deltad_2D)),(abs.(deltav_2D)), bins=(10 .^(-3:0.2:10), 10 .^(-3:.2:10)), cmap="afmhot_r")
#plot((x), (y))
#cbar=colorbar(orientation="vertical")
#ylabel(L"(\delta v / <v>)_{2D}", fontsize=15)
#xlabel(L"(\delta \rho / <\rho>)_{2D}", fontsize=15)
#xticks(fontsize=15)
#yticks(fontsize=15)
#yscale("log")
#xscale("log")
#xlim(-3,0)
#ylim(-3,0)

#readline()
#clf()

println("-----------  2D FLUCTUATIONS: temperature-velocity -------------")

m(x, p) =  p[1] .* x


p0=[1.]
fit = curve_fit(m, (abs.(deltaT_2D)), (abs.(deltav_2D)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid

x = collect(minimum(abs.(deltaT_2D)):0.001:maximum(abs.(deltaT_2D)))
y = zeros(length(x))

for i in 1:length(x)
    #y[i] = param[1] + param[2]*log10(abs(x[i]))
    y[i] = param[1]*((x[i]))
end

R = Statistics.cor(abs.(deltaT_2D), abs.(deltav_2D))
p_value=pvalue(CorrelationTest(abs.(deltaT_2D),abs.(deltav_2D)))

#ans = "2D FLUCTUATIONS temperature-velocity:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/02_radial_filtering/new_",name[l],"_2D_radial_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1=string("/home/marco/Scrivania/Ettori_project/02_2D_radial_filtering/x_2D_02_radial_filtering_temperature_velocity_fit_parameters.txt")
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end


#hist, xedges, yedges=hist2D((abs.(deltaT_2D)),(abs.(deltav_2D)), bins=(10 .^(-3:0.2:0), 10 .^(-3:.2:0)), cmap="afmhot_r")
#plot((x), (y))
#cbar=colorbar(orientation="vertical")
#ylabel(L"(\delta v / <v>)_{2D}", fontsize=15)
#xlabel(L"(\delta T / <T>)_{2D}", fontsize=15)
#xticks(fontsize=15)
#yticks(fontsize=15)
#yscale("log")
#xscale("log")
#xlim(-3,0)
#ylim(-3,0)

#readline()
#clf()

println("----------- 3D FLUCTUATIONS: temperature-velocity^2 -------------")

m(x, p) = p[1] .* x
p0=[1.]
fit = curve_fit(m, (abs.(deltaT_2D)), (abs.(deltav_2D)).^2, p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid

R = Statistics.cor(abs.(deltaT_2D), deltav_2D.^2)
p_value=pvalue(CorrelationTest(abs.(deltaT_2D),deltav_2D.^2))

#ans = "2D FLUCTUATIONS temperature-velocity^2:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/02_radial_filtering/new_",name[l],"_2D_radial_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1="/home/marco/Scrivania/Ettori_project/02_2D_radial_filtering/x_2D_02_radial_temperature_velocity2_filtering_fit_parameters.txt"
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end


x = collect(minimum(abs.(deltaT_2D)):0.001:maximum(abs.(deltaT_2D)))
y = zeros(length(x))

for i in 1:length(x)
    y[i] = param[1]*(abs(x[i]))
    #y[i] = param[1]*(abs(x[i]))
end

#hist, xedges, yedges=hist2D((abs.(deltaT_2D)),(abs.(deltav_2D)).^2, bins=(10 .^(-3:0.2:0.5), 10 .^(-3:.2:0.5)), cmap="afmhot_r")
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


end
