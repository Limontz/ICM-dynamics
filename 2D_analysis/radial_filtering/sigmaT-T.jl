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

main = "/home/marco/Scrivania/Ettori_project/Codes/2D_analysis/radial_filtering/"
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
include(string(main, "y_Projection.jl"))
include(string(main, "vel_module2D.jl"))
include(string(main, "2D_radial_profile.jl"))

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

    dx=0.02
    d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D = projection(d, temp, vx, vy, vz, Npoint)
    #v_2D = vmodule2D(vx,vy,vz,Npoint)

    #********************* Let's consider only ************************
    #********************* the cells inside R500 **********************
    Ncell = floor(Int32, (r500/(20*kpc)))

    N = 15
    Ncell = floor(Int32, (r500/(20*kpc)) + N)


    d_2D = d_2D[xc-Ncell-N:xc+Ncell+N,
                yc-Ncell-N:yc+Ncell+N]

    temp_2D = spec_T_2D[xc-Ncell-N:xc+Ncell+N,
                        yc-Ncell-N:yc+Ncell+N]

    vx_2D = vx_2D[xc-Ncell-N:xc+Ncell+N,
                  yc-Ncell-N:yc+Ncell+N]

    vy_2D = vy_2D[xc-Ncell-N:xc+Ncell+N,
               yc-Ncell-N:yc+Ncell+N]

    vz_2D = vz_2D[xc-Ncell-N:xc+Ncell+N,
               yc-Ncell-N:yc+Ncell+N]

    d_max, index = findmax(d_2D)
    xc = index[1]
    yc = index[2]

    xc = floor(Int32, length(d_2D[:,1])/2 +1)
    yc = floor(Int32, length(d_2D[1,:])/2 +1)

    Npoint=length(d_2D[1,:])
    Nhalf = Int(floor(Npoint/2))
    r_max = 20. * kpc * Ncell
    #r_max=r500
    rmin = 20. * kpc
    radial_bin = 10

    d_2D_prof, temp_2D_prof, vx_2D_prof, vy_2D_prof, vz_2D_prof, radius = profile2D(d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D, xc, yc, zc, r500, Nhalf, radial_bin, thickness)

    Ncell = floor(Int32, (r500/(20*kpc)))

    d_2D = d_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
    spec_T_2D = spec_T_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
    vx_2D = vx_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
    vy_2D = vy_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
    vz_2D = vz_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]

    sigma_T[l] = std(spec_T_2D)*8.6173e-8
    sigma_d[l] = std(d_2D)
    #sigma_v[l] = std(v)

    sigma_T_prof[l] = abs(temp_2D_prof[length(snap)] - temp_2D_prof[1])*8.6173e-8
    sigma_d_prof[l] = abs(d_2D_prof[length(snap)] - d_2D_prof[1])
    #sigma_v_prof[l] = v_prof[length(snap)] - v_prof[1]

    T_mean[l] = mean(spec_T_2D)*8.6173e-8

end

output1=file1=string("/home/marco/Scrivania/Ettori_project/Plots/y_sigmaT-sigman_2D.txt")
out1=open(output1,"w") do out1
writedlm( out1, [sigma_T sigma_d T_mean])
end
