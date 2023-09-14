using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit


main = "/home/marco/Scrivania/Ettori_project/Codes/2D_analysis/radial_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"
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
dr = 0.1*r500 #thick of shell for radial profile
#r_max = 160 * 20 *kpc * sqrt(3)
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

Ncell = floor(Int32, (r500/(20*kpc))) + 10

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

include(string(main, "Projection.jl"))
include(string(main, "vel_module2D.jl"))
dx=0.02
d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D = projection(d, temp, vx, vy, vz, Npoint)
#v_2D = vmodule2D(vx,vy,vz,Npoint)
pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((d_2D))))

axis("scaled")
#title("IT90_3", fontsize= "20")
xticks(fontsize= 20)
yticks(fontsize=20);
cbar=colorbar(orientation="vertical")
cbar.ax.tick_params(labelsize=20)
cbar.ax.yaxis.get_offset_text().set(size=20)
#cbar[:set_yticklabel](labelsize="large")
#cbar[:set_label](L"Log(T [k])",fontsize= 25)
xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 30)
ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 30)
readline()
clf()
#****************** RADIAL PROFILE **println("r=",r)**************************
include(string(main, "radial_profile.jl"))

Npoint=length(d_2D[1,:])
Nhalf = Int(floor(Npoint/2))
r_max = 20. * kpc * Ncell * sqrt(3)
#r_max=r500
rmin = 20. * kpc
radial_bin = Int(floor(10*r_max/r500))



include(string(main, "2D_radial_profile.jl"))
d_2D_prof, temp_2D_prof, vx_2D_prof, vy_2D_prof, vz_2D_prof, radius = profile2D(d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D, xc, yc, zc, r500, Nhalf, radial_bin)

#****************** FILTERING *********************************

include(string(main, "2D_radial_filtering.jl"))
include(string(main, "v_mean2D.jl"))
deltad_2D, deltaT_2D, deltav_2D = filtering2D(d_2D, d_2D_prof, spec_T_2D, temp_2D_prof,
                                              vx_2D, vx_2D_prof, vy_2D, vy_2D_prof, vz_2D, vz_2D_prof,
                                              radial_bin, Npoint, xc,yc,zc, r500)


println("2D filtering done")

Ncell = floor(Int32, (r500/(20*kpc)))

deltad_2D = deltad_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
deltaT_2D = deltaT_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
deltav_2D = deltav_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]

Npoint = length(deltad_2D[1,:])

dx=0.02
subplot(1,2,1)
pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((deltad_2D))))

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

subplot(1,2,2)
pcolormesh((-Ncell*dx:dx:Ncell*dx),(-Ncell*dx:dx:Ncell*dx),(((deltaT_2D))))

axis("scaled")
#title("IT90_3", fontsize= "20")
xticks(fontsize= 20)
yticks(fontsize=20);
cbar=colorbar(orientation="vertical")
cbar.ax.tick_params(labelsize=20)
cbar.ax.yaxis.get_offset_text().set(size=20)
#cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"Log(\delta T)",fontsize= 25)
xlabel(string("y ",L"[\mathrm{Mpc}]"),fontsize= 30)
ylabel(string("z ",L"[\mathrm{Mpc}]"),fontsize= 30)

readline()
clf()

#************************ masking the largest fluctuations ***********************

include(string(main, "clump_masking_2D.jl"))
include(string(main, "pdf.jl"))
include(string(main, "large_fluctuations_mask_2D.jl"))


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

#
#plot(deltad_2D, deltaT_2D, linestyle="", marker=".")
#hjj
#output1=string("/home/marco/Scrivania/Ettori_project/test.txt")
#out1=open(output1,"a+") do out1
#writedlm( out1, [deltaT_2D, deltad_2D, deltav_2D])
#end

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


println("----------- 2D FLUCTUATIONS: temperature-density -------------")
println("m = ", param[1], " +- ", sigma_par[1])
#println(L"m = ", param[2], " +- ", sigma_par[2])

x = collect(minimum((deltad_2D)):0.001:maximum((deltad_2D)))
y = zeros(length(x))

for i in 1:length(x)
    #y[i] = param[1] + param[2]*log10(abs(x[i]))
    y[i] = param[1]*((x[i]))
end

deltad2D_mean = mean((deltad_2D))
deltaT2D_mean = mean((deltaT_2D))
somma = 0.
num = 0
den = 0
for i in 1:length(deltad_2D)
    global somma = somma + ((deltaT_2D[i]) - param[1] * (deltad_2D[i]))^2 / (param[1] * (deltad_2D[i]))
    global num = num + ((deltaT_2D[i]) - param[1]*(deltad_2D[i]) )^2
    global den = den + ((deltaT_2D[i]) - deltaT2D_mean)^2
end


num1 = ((deltad_2D[1]) - mean((deltad_2D))) * ((deltaT_2D[1]) - deltaT2D_mean)
den1 = ((deltad_2D[1]) - mean((deltad_2D)))^2
den2 = ((deltaT_2D[1]) - deltaT2D_mean)^2

for i in 2:length(deltad_2D)

    global num1 = num1 + ((deltad_2D[i]) - deltad2D_mean) * ((deltaT_2D[i]) - deltaT2D_mean)
    global den1 = den1 + ((deltad_2D[i]) - deltad2D_mean)^2
    global den2 = den2 + ((deltaT_2D[i]) - deltaT2D_mean)^2

end
R = num1 / ( sqrt(den1 * den2))
println("R squared =", 1 - num/den)
#println("Pearson coefficient =", sqrt(1 - num/den))
println("Pearson coefficient =", R)
println("-----------------------------------------")

#h, xedges, yedges=hist2D(log10.(abs.(deltad_2D)),log10.(abs.(deltaT_2D)), bins=(30,30), cmap="afmhot_r")
#hist, xedges, yedges=hist2D(((deltad_2D)),((deltaT_2D)), bins=(30,30), cmap="afmhot_r")
#hist=hist2D((abs.(deltad_2D)),(abs.(deltaT_2D)), bins=(xedges,yedges), cmap="afmhot_r")
plot(deltad_2D, deltaT_2D, linestyle="", marker=".")
plot((x), (y))

#cbar=colorbar(orientation="vertical")
ylabel(L"(\delta T /<T>)_{2D}", fontsize=15)
xlabel(L"(\delta \rho / <\rho>)_{2D}", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)
#yscale("log")
#xscale("log")
#xlim(-0.3,0.3)
#ylim(-0.3,0.3)
readline()
clf()


println("-----------  2D FLUCTUATIONS: density-velocity -------------")

m(x, p) =  p[1] .+ p[2].* x
p0 = [0.,1.]

fit = curve_fit(m, (abs.(deltad_2D)), (abs.(deltav_2D)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid

x = collect(minimum(abs.(deltad_2D)):0.001:maximum(abs.(deltad_2D)))
y = zeros(length(x))

for i in 1:length(x)
    #y[i] = param[1] + param[2]*log10(abs(x[i]))
    y[i] = param[1] + param[2]*(abs(x[i]))
end

deltav2D_mean = mean(abs.(deltav_2D))
somma = 0.
num = 0
den = 0
for i in 1:length(deltad_2D)
    #global somma = somma + (abs.(deltav_2D[i]) - param[1] * abs.(deltad_2D[i]))^2 / (param[1] * abs.(deltad_2D[i]))
    global num = num + (abs.(deltav_2D[i]) - param[1]*abs.(deltad_2D[i]) )^2
    global den = den + (abs.(deltav_2D[i]) - deltav2D_mean)^2
end

num1 = zeros(length(deltad_2D) )
den1 = zeros(length(deltad_2D))
den2 = zeros(length(deltad_2D))
num1[1] = (abs.(deltad_2D[1]) - mean(abs.(deltad_2D))) * (abs.(deltav_2D[1]) - mean(abs.(deltav_2D)))
den1[1] = (abs.(deltad_2D[1]) - mean(abs.(deltad_2D)))^2
den2[1] = (abs.(deltav_2D[1]) - mean(abs.(deltav_2D)))^2

for i in 2:length(deltad_2D)

    num1[i] = num1[i-1] + (abs.(deltad_2D[i]) - deltad2D_mean) * (abs.(deltav_2D[i]) - deltav2D_mean)
    den1[i] = den1[i-1] + (abs.(deltad_2D[i]) - deltad2D_mean)^2
    den2[i] = den2[i-1] + (abs.(deltav_2D[i]) - deltav2D_mean)^2

end
R = num1[length(num1)] / ( sqrt(den1[length(den1)] * den2[length(den2)]))

println("m = ", param[2], " +- ", sigma_par[2])
#println(L"m = ", param[2], " +- ", sigma_par[2])
println("R squared =", 1 - num/den)
#println("Pearson coefficient =", sqrt(1 - num/den))
println("Pearson coefficient =", R)
println("-----------------------------------------")


hist, xedges, yedges=hist2D((abs.(deltad_2D)),(abs.(deltav_2D)), bins=(10 .^(-3:0.2:10), 10 .^(-3:.2:10)), cmap="afmhot_r")
plot((x), (y))
cbar=colorbar(orientation="vertical")
ylabel(L"(\delta v / <v>)_{2D}", fontsize=15)
xlabel(L"(\delta \rho / <\rho>)_{2D}", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)
yscale("log")
xscale("log")
#xlim(-3,0)
#ylim(-3,0)

readline()
clf()

println("-----------  2D FLUCTUATIONS: temperature-velocity -------------")

m(x, p) =  p[1] .+ p[2].* x


p0=[0.,1.]
fit = curve_fit(m, (abs.(deltaT_2D)), (abs.(deltav_2D)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid

x = collect(minimum(abs.(deltaT_2D)):0.001:maximum(abs.(deltaT_2D)))
y = zeros(length(x))

for i in 1:length(x)
    #y[i] = param[1] + param[2]*log10(abs(x[i]))
    y[i] = param[1] + param[2]*((x[i]))
end

deltav2D_mean = mean(abs.(deltav_2D))
deltaT2D_mean = mean(abs.(deltaT_2D))
somma = 0.
num = 0
den = 0
for i in 1:length(deltaT_2D)
    #global somma = somma + (abs.(deltav_2D[i]) - param[1] * abs.(deltad_2D[i]))^2 / (param[1] * abs.(deltad_2D[i]))
    global num = num + (abs.(deltav_2D[i]) - param[1]*abs.(deltaT_2D[i]) )^2
    global den = den + (abs.(deltav_2D[i]) - deltav2D_mean)^2
end

num1 = zeros(length(deltaT_2D) )
den1 = zeros(length(deltaT_2D))
den2 = zeros(length(deltaT_2D))
num1[1] = (abs.(deltaT_2D[1]) - mean(abs.(deltaT_2D))) * (abs.(deltav_2D[1]) - mean(abs.(deltav_2D)))
den1[1] = (abs.(deltaT_2D[1]) - mean(abs.(deltaT_2D)))^2
den2[1] = (abs.(deltav_2D[1]) - mean(abs.(deltav_2D)))^2

for i in 2:length(deltaT_2D)

    num1[i] = num1[i-1] + (abs.(deltaT_2D[i]) - deltaT2D_mean) * (abs.(deltav_2D[i]) - deltav2D_mean)
    den1[i] = den1[i-1] + (abs.(deltaT_2D[i]) - deltaT2D_mean)^2
    den2[i] = den2[i-1] + (abs.(deltav_2D[i]) - deltav2D_mean)^2

end
R = num1[length(num1)] / ( sqrt(den1[length(den1)] * den2[length(den2)]))

println("m = ", param[2], " +- ", sigma_par[2])
#println(L"m = ", param[2], " +- ", sigma_par[2])
println("R squared =", 1 - num/den)
#println("Pearson coefficient =", sqrt(1 - num/den))
println("Pearson coefficient =", R)
println("-----------------------------------------")


hist, xedges, yedges=hist2D((abs.(deltaT_2D)),(abs.(deltav_2D)), bins=(10 .^(-3:0.2:0), 10 .^(-3:.2:0)), cmap="afmhot_r")
plot((x), (y))
cbar=colorbar(orientation="vertical")
ylabel(L"(\delta v / <v>)_{2D}", fontsize=15)
xlabel(L"(\delta T / <T>)_{2D}", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)
yscale("log")
xscale("log")
#xlim(-3,0)
#ylim(-3,0)

readline()
clf()

fx






using KernelDensity
using PyCall
using PyPlot

#dens = kde((deltaT,deltav))
#println(typeof(dens))

bins = 5
deltaT1=[]
deltav1=[]
deltad1=[]
deltaT2=[]
deltav2=[]
deltad2=[]
deltaT3=[]
deltav3=[]
deltad3=[]
deltaT4=[]
deltav4=[]
deltad4=[]
deltaT5=[]
deltav5=[]
deltad5=[]
deltaT6=[]
deltav6=[]
deltad6=[]

volume = Ncell^3/5
println(((volume*5)^(1/3)* 20 *kpc)/r500)
rad = zeros(bins+1)


for i in 1:bins+1

    rad[i] = ((i-1)*volume)^(1/3)

end
for i in 1:bins
    println(rad[i+1]^3-rad[i]^3)
end
#println(rad)
#println((rad * 20 *kpc)/r500)

@inbounds for i in xc - Ncell: xc + Ncell #(xc-N_half:xc+N_half) se voglio una boxcon Ncell centrata in xc,yc,zc)
@inbounds @simd for j in yc - Ncell: yc + Ncell
  #@inbounds @simd for k in zc - Ncell: zc + Ncell

            @fastmath r=Int(floor(sqrt((j-yc -0.5)^2. +(i- xc -0.5)^2.)))
                      ra=r#convert(Int64,trunc(r))  # assign radial bin
                      if ra <= minimum(rad)
                         ra=minimum(rad)
                       end
                       if ra >= maximum(rad)
                          ra=maximum(rad)
                       end

                       #println(ra, " ", rad)

                       if (ra > rad[1]-1 && ra <= rad[2]-1)
                          append!(deltaT1, deltaT[i,j,k])
                          append!(deltav1, deltav[i,j,k])
                          append!(deltad1, deltad[i,j,k])
                       end

                       if (ra > rad[2]-1 && ra <= rad[3]-1)
                          append!(deltaT2, deltaT[i,j,k])
                          append!(deltav2, deltav[i,j,k])
                          append!(deltad2, deltad[i,j,k])
                       end

                       if (ra > rad[3]-1 && ra <= rad[4]-1)
                          append!(deltaT3, deltaT[i,j,k])
                          append!(deltav3, deltav[i,j,k])
                          append!(deltad3, deltad[i,j,k])
                       end

                       if (ra > rad[4]-1 && ra <= rad[5]-1)
                          append!(deltaT4, deltaT[i,j,k])
                          append!(deltav4, deltav[i,j,k])
                          append!(deltad4, deltad[i,j,k])
                       end

                       if (ra > rad[5]-1 && ra <= rad[6]-1)
                          append!(deltaT5, deltaT[i,j,k])
                          append!(deltav5, deltav[i,j,k])
                          append!(deltad5, deltad[i,j,k])
                       end

                  end
              end
          end
#print(length(deltad1), " ",length(deltaT2)," ",length(deltaT3)," ",length(deltaT4)," ",length(deltaT5))
#plot(abs.(deltad1),abs.(deltav1), linestyle="", marker=".", label=L"0-0.5 \ R_{500}")
#plot(abs.(deltad2),abs.(deltav2), linestyle="", marker=".", label=L"0.5-0.7 \ R_{500}")
#plot(abs.(deltad3),abs.(deltav3), linestyle="", marker=".", label=L"0.7-0.8 \ R_{500}")
#plot(abs.(deltad4),abs.(deltav4), linestyle="", marker=".", label=L"0.8-0.9 \ R_{500}")
#plot(abs.(deltad5),abs.(deltav5), linestyle="", marker=".", label=L"0.9-1 \ R_{500}")
#legend()

#pcolormesh(abs.(dens.x),dens.y,log10.(dens.density), cmap="binary")
xlabel(L"\delta \rho / <\rho>", fontsize=15)
ylabel(L"\mathcal{M}", fontsize=15)
xticks(fontsize=15)
yticks(fontsize=15)
yscale("log")
xscale("log")
