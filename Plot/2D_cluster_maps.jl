using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit

main = "/home/marco/Scrivania/Ettori_project/Codes/2D_analysis/box_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"
snap = "0187"

const mu = 0.62
const mp = 1.67e-24
const kb = 1.28e-16
const pc = 3.08e18
const kpc = (10.0^3) * pc
const G = 6.67e-8
const gamma = 5.0 / 3.0

rvir=0.88*1.e3*kpc
r500 = 5*rvir/7
r_max=r500

include(string(main, "read.jl"))

d, d_dm, temp, vx, vy, vz, Npoint, conv_d  = read(snap, cluster)
#v = vmodule3D(vx,vy,vz,Npoint)

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

Ncell = floor(Int32, (r500/(20*kpc)))
include(string(main, "z_Projection.jl"))
dx=0.02

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

N = 0 # floor(Int32, thickness*200)
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

d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D = projection(d, temp, vx, vy, vz, Npoint)

Nhalf_box = 0
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

v_2D = zeros(length(vx_2D[:,1]), length(vx_2D[:,1]))
for i in 1:length(vx_2D[:,1])
    for j in 1:length(vx_2D[:,1])

        v_2D[i,j] = sqrt(vx_2D[i,j]^2 + vy_2D[i,j]^2 + vz_2D[i,j]^2)

    end
end

xc = floor(Int32, length(d_2D[:,1])/2 +1)
yc = floor(Int32, length(d_2D[1,:])/2 +1)



# Create the figure and subplots
fig, axes = subplots(1, 3, figsize=(15, 5))

# First subplot
im1 = axes[1].pcolormesh((-(Ncell + Nhalf_box) * dx:dx:(Ncell + Nhalf_box) * dx), (-(Ncell + Nhalf_box) * dx:dx:(Ncell + Nhalf_box) * dx), log10.(d_2D))
axes[1].axis("scaled")
axes[1].set_xlabel("y [Mpc]", fontsize=16)
axes[1].set_ylabel("z [Mpc]", fontsize=16)
axes[1].tick_params(axis="x", labelsize=16)
axes[1].tick_params(axis="y", labelsize=16)
cbar1 = colorbar(im1, ax=axes[1], fraction=0.05, pad=0.01)  # Adjust height as needed
cbar1.ax.tick_params(labelsize=13)
cbar1.ax.yaxis.get_offset_text().set(size=13)
cbar1.set_label(L"Log (\rho \ [g/cm^3])", fontsize=16)

# Second subplot
im2 = axes[2].pcolormesh((-(Ncell + Nhalf_box) * dx:dx:(Ncell + Nhalf_box) * dx), (-(Ncell + Nhalf_box) * dx:dx:(Ncell + Nhalf_box) * dx), temp_2D*8.62e-8)
axes[2].axis("scaled")
axes[2].set_xlabel("y [Mpc]", fontsize=16)
axes[2].set_ylabel("z [Mpc]", fontsize=16)
cbar2 = colorbar(im2, ax=axes[2], fraction=0.05, pad=0.01) # Adjust height as needed
cbar2.ax.tick_params(labelsize=13)
cbar2.ax.yaxis.get_offset_text().set(size=13)
cbar2.set_label(L"T [keV]", fontsize=16)
axes[2].tick_params(axis="x", labelsize=16)
axes[2].tick_params(axis="y", labelsize=16)

# Third subplot
im3 = axes[3].pcolormesh((-(Ncell + Nhalf_box) * dx:dx:(Ncell + Nhalf_box) * dx), (-(Ncell + Nhalf_box) * dx:dx:(Ncell + Nhalf_box) * dx), log10.(v_2D*1.e-5))
axes[3].axis("scaled")
axes[3].set_xlabel("y [Mpc]", fontsize=16)
axes[3].set_ylabel("z [Mpc]", fontsize=16)
cbar3 = colorbar(im3, fraction=0.05, pad=0.01)  # Adjust height as needed
cbar3.ax.tick_params(labelsize=13)
cbar3.ax.yaxis.get_offset_text().set(size=13)
cbar3.set_label(L"Log(v [km/s])", fontsize=16)
axes[3].tick_params(axis="x", labelsize=16)
axes[3].tick_params(axis="y", labelsize=16)

tight_layout()

show()
