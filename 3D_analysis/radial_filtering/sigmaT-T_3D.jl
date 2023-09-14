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

main = "/home/marco/Scrivania/Ettori_project/Codes/3D_analysis/radial_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"

const mu = 0.62
const mp = 1.67e-24
const kb = 1.28e-16
const pc = 3.08e18
const kpc = (10.0^3) * pc
const G = 6.67e-8
const gamma = 5.0 / 3.0
thickness = 0.1

include(string(main, "read.jl"))
include(string(main, "vel_module3D.jl"))
include(string(main, "nfluct.jl"))
include(string(main, "radial_profile.jl"))

sigma_T = zeros(length(snap))
sigma_d = zeros(length(snap))
#sigma_v = zeros(length(snap))

sigma_T_prof = zeros(length(snap))
sigma_d_prof = zeros(length(snap))
#sigma_v_prof = zeros(length(snap))

T_mean = zeros(length(snap))

for l in 1:length(snap)

    rvir=r_vir[l]*1.e3*kpc
    r500 = 5*rvir/7

    d, d_dm, temp, vx, vy, vz, Npoint, conv_d  = read(snap[l], cluster)

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

    Npoint=length(d[1,1,:])
    Nhalf = Int(floor(Npoint/2))

    radial_bin = 10
    rmin = 20. * kpc

    d_prof, temp_prof,vx_prof, vy_prof, vz_prof, radius =profile(d, temp, vx, vy, vz, xc, yc, zc, r500, rmin, Nhalf,radial_bin, thickness)##

    Ncell = floor(Int32, (r500/(20*kpc)))

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

    sigma_T[l] = std(temp)*8.6173e-8
    sigma_d[l] = std(d)
    #sigma_v[l] = std(v)

    sigma_T_prof[l] = abs(temp_prof[length(snap)] - temp_prof[1])*8.6173e-8
    sigma_d_prof[l] = abs(d_prof[length(snap)] - d_prof[1])
    #sigma_v_prof[l] = v_prof[length(snap)] - v_prof[1]

    T_mean[l] = mean(temp)*8.6173e-8

    #println(sigma_T[l], " ", sigma_d[l], " ",T_mean[l])


end

output1=file1=string("/home/marco/Scrivania/Ettori_project/Plots/sigmaT-sigman.txt")
out1=open(output1,"w") do out1
writedlm( out1, [sigma_T sigma_d T_mean])
end
