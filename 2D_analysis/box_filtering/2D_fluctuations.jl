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

main = "/home/marco/Scrivania/Ettori_project/Codes/2D_analysis/box_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"

include(string(main, "2D_box_filtering.jl"))
include(string(main, "v_mean2D.jl"))
include(string(main, "z_Projection.jl"))
include(string(main, "vel_module2D.jl"))
include(string(main, "read.jl"))
include(string(main, "vel_module3D.jl"))
include(string(main, "nfluct.jl"))

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

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

N = 15
Ncell = floor(Int32, (r500/(20*kpc)) + N)

d_dm = d_dm[xc-Ncell-N:xc+Ncell+N,
            yc-Ncell-N:yc+Ncell+N,
            zc-Ncell-N:zc+Ncell+N]


d = d[xc-Ncell-N:xc+Ncell+N,
            yc-Ncell-N:yc+Ncell+N,
            zc-Ncell-N:zc+Ncell+N]

temp = temp[xc-Ncell-N:xc+Ncell+N,
            yc-Ncell-N:yc+Ncell+N,
            zc-Ncell-N:zc+Ncell+N]

vx = vx[xc-Ncell-N:xc+Ncell+N,
            yc-Ncell-N:yc+Ncell+N,
            zc-Ncell-N:zc+Ncell+N]

vy = vy[xc-Ncell-N:xc+Ncell+N,
            yc-Ncell-N:yc+Ncell+N,
            zc-Ncell-N:zc+Ncell+N]

vz = vz[xc-Ncell-N:xc+Ncell+N,
            yc-Ncell-N:yc+Ncell+N,
            zc-Ncell-N:zc+Ncell+N]

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

Npoint=length(d[:,1,1])

#*********************PROJECTION***********************************

Ncell = floor(Int32, (r500/(20*kpc)))

dx=0.02
d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D = projection(d, temp, vx, vy, vz, Npoint)
#v_2D = vmodule2D(vx,vy,vz,Npoint)

#********************* Let's consider only ************************
#********************* the cells inside R500 **********************
Ncell = floor(Int32, (r500/(20*kpc)))
Ncell_box= 15
Nhalf_box=floor(Int32,Ncell_box/2)


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

d_max, index = findmax(d_2D)
xc = index[1]
yc = index[2]

xc = floor(Int32, length(d_2D[:,1])/2 +1)
yc = floor(Int32, length(d_2D[1,:])/2 +1)


#dx=0.02
#subplot(1,2,1)
#pcolormesh(((-Ncell-Nhalf_box)*dx:dx:(Ncell+Nhalf_box)*dx),((-Ncell-Nhalf_box)*dx:dx:(Ncell+Nhalf_box)*dx),(((vx_2D))))

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


deltad_2D, deltaT_2D, deltav_2D = box_filtering2D(d_2D, temp_2D, vx_2D, vy_2D, vz_2D, Ncell+Nhalf_box, Nhalf_box)
println("2D filtering done")

deltad_2D = deltad_2D[xc-Ncell:xc+Ncell,
            yc-Ncell:yc+Ncell]

deltaT_2D = deltaT_2D[xc-Ncell:xc+Ncell,
                    yc-Ncell:yc+Ncell]

deltav_2D = deltav_2D[xc-Ncell:xc+Ncell,
              yc-Ncell:yc+Ncell]

d_max, index = findmax(d_2D)
xc = index[1]
yc = index[2]

sigmaT[l] = std(deltaT_2D)
sigman[l] = std(deltad_2D)


#println("length d_2d= ", length(d_2D[1,:]) )

#println("xc based on d2d= ", xc)

#xc = floor(Int32, length(d_2D[:,1])/2 +1)
#println("xc based on len/2= ", xc)
#yc = floor(Int32, length(d_2D[1,:])/2 +1)

#println("2D filtering done")

#subplot(1,2,2)
#pcolormesh(((-Ncell)*dx:dx:(Ncell)*dx),((-Ncell)*dx:dx:(Ncell)*dx),(((deltaT_2D))))

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

#****************** DeltaT/T-deltav relation *********************************

Ncell = floor(Int32, (r500/(20*kpc)))
println((Ncell), " ", 0.6*r500/kpc)

#-------------- Fluctuations as result of filtering --------------------------------

deltaT_2D = vec(deltaT_2D)
deltad_2D = vec(deltad_2D)
deltav_2D = vec(deltav_2D)

println("vectorization done")
#************************* LINEAR FIT ********************************************
println("----------- 2D FLUCTUATIONS: temperature-density -------------")


m(x, p) =  p[1] .* x
#m(x, p) = log.(abs.(1 .+ p[1] .* x[:,2].^4 .+ p[2] .* x[:,2].^2 .* x[:,1] .* x[:,2]) )
#m(x,p) = log.(1. .+ p[1]^2 .* x.^2)
p0=[1.]
fit = curve_fit(m, ((deltad_2D)), ((deltaT_2D)), p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid


println("----------- 2D FLUCTUATIONS: temperature-density -------------")


x = collect(minimum((deltad_2D)):0.001:maximum((deltad_2D)))
y = zeros(length(x))

for i in 1:length(x)
    #y[i] = param[1] + param[2]*log10(abs(x[i]))
    y[i] = param[1]*((x[i]))
end

#h, xedges, yedges=hist2D(log10.(abs.(deltad_2D)),log10.(abs.(deltaT_2D)), bins=(30,30), cmap="afmhot_r")
#plot(deltad_2D, deltaT_2D, linestyle="", marker=".")
#hist=hist2D((abs.(deltad_2D)),(abs.(deltaT_2D)), bins=(xedges,yedges), cmap="afmhot_r")
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
#readline()
#clf()

R = Statistics.cor(deltad_2D, deltaT_2D)
p_value=pvalue(CorrelationTest(deltad_2D,deltaT_2D))

#ans = "2D FLUCTUATIONS density-temperature:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/100_kpc_filtering/new_",name[l],"_2D_box_filtering_fit_parameters.txt")
#out1=open(output1,"w") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1=string("/home/marco/Scrivania/Ettori_project/300_kpc_2D_box_filtering/z_2D_300_kpc_box_filtering_density_temperature_fit_parameters.txt")
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end

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


#hist, xedges, yedges=hist2D((abs.(deltad_2D)),(abs.(deltav_2D)), bins=(10 .^(-3:0.2:1), 10 .^(-3:.2:1)), cmap="afmhot_r")
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

R = Statistics.cor(abs.(deltad_2D), abs.(deltav_2D))
p_value=pvalue(CorrelationTest(abs.(deltad_2D),abs.(deltav_2D)))

#ans = "2D FLUCTUATIONS density-velocity:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/100_kpc_filtering/new_",name[l],"_2D_box_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1=string("/home/marco/Scrivania/Ettori_project/300_kpc_2D_box_filtering/z_2D_300_kpc_box_filtering_density_velocity_fit_parameters.txt")
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end

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

R = Statistics.cor(abs.(deltaT_2D), abs.(deltav_2D))
p_value=pvalue(CorrelationTest(abs.(deltaT_2D), abs.(deltav_2D)))

#ans = "2D FLUCTUATIONS temperature-velocity:"
#output1=string("/home/marco/Scrivania/Ettori_project/100_kpc_2D_box_filtering/x_2D_100_kpc_box_filtering_temperature_velocity_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1=string("/home/marco/Scrivania/Ettori_project/300_kpc_2D_box_filtering/z_2D_300_kpc_box_filtering_temperature_velocity_fit_parameters.txt")
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end

println("----------- 3D FLUCTUATIONS: temperature-velocity^2 -------------")

p0=[1.]
fit = curve_fit(m, (abs.(deltaT_2D)), (abs.(deltav_2D)).^2, p0)
param= fit.param
sigma_par= stderror(fit)
res= fit.resid

println(param)
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
#xlabel(L"(\delta T / <T>)_{
#xticks(fontsize=15)
#yticks(fontsize=15)
#yscale("log")
#xscale("log")
#xlim(-3,0)
#ylim(-3,0)
#readline()

R = Statistics.cor(abs.(deltaT_2D), deltav_2D.^2)
p_value=pvalue(CorrelationTest(abs.(deltaT_2D),deltav_2D.^2))

#ans = "2D FLUCTUATIONS temperature-velocity^2:"
#output1=string("/home/marco/Scrivania/Ettori_project/",name[l],"/100_kpc_filtering/new_",name[l],"_2D_box_filtering_fit_parameters.txt")
#out1=open(output1,"a+") do out1
#write( out1, ans)
#end
#close(output1)

output1=file1=string("/home/marco/Scrivania/Ettori_project/300_kpc_2D_box_filtering/z_2D_300_kpc_box_filtering_temperature_velocity2_fit_parameters.txt")
out1=open(output1,"a+") do out1
writedlm( out1, [param[1] sigma_par[1] R p_value])
end

end
