using HDF5
using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using LsqFit
@pyimport mpl_toolkits.axes_grid1 as axgrid


main = "/home/marco/Scrivania/Ettori_project/Codes/2D_analysis/radial_filtering/"
cluster =  "/home/marco/Scrivania/Tesi/Ammassi_tesi/"
snap = "0187"

const mu = 0.62
const mp = 1.67e-24
const kb = 1.28e-16
const pc = 3.08e18
const kpc = (10.0^3) * pc
const G = 6.67e-8
const gamma = 5.0 / 3.0

include(string(main, "read.jl"))
include(string(main, "z_Projection.jl"))
include(string(main, "2D_radial_profile.jl"))
include(string(main, "2D_radial_filtering.jl"))
include(string(main, "v_mean2D.jl"))
include(string(main, "../box_filtering/2D_box_filtering.jl"))
include(string(main, "../box_filtering/v_mean2D.jl"))

rvir=0.88*1.e3*kpc
r500 = 5*rvir/7
r_max=r500
thickness = 0.1


d, d_dm, temp, vx, vy, vz, Npoint, conv_d  = read(snap, cluster)
#v = vmodule3D(vx,vy,vz,Npoint)

d_max, index = findmax(d_dm)
xc = index[1]
yc = index[2]
zc = index[3]

Ncell = floor(Int32, (r500/(20*kpc)))
dx=0.02

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


d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D = projection(d, temp, vx, vy, vz, Npoint)

v_2D = zeros(length(vx_2D[:,1]),length(vx_2D[:,1]))
for i in 1:length(vx_2D[:,1])
    for j in 1:length(vx_2D[1,:])

        v_2D[i,j] = sqrt(vx_2D[i,j]^2 + vy_2D[i,j]^2 + vz_2D[i,j]^2)

    end
end


Npoint=length(d_2D[1,:])
Nhalf = Int(floor(Npoint/2))
r_max = 20. * kpc * Ncell
radial_bin = Int(floor((1/thickness)*r_max/r500 + 1))


d_2D_prof, temp_2D_prof, vx_2D_prof, vy_2D_prof, vz_2D_prof, radius = profile2D(d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D, xc, yc, zc, r500, Nhalf, radial_bin, thickness)

radial_deltad_2D, radial_deltaT_2D, radial_deltav_2D = filtering2D(d_2D, d_2D_prof, spec_T_2D, temp_2D_prof,
                                              vx_2D, vx_2D_prof, vy_2D, vy_2D_prof, vz_2D, vz_2D_prof,
                                              radial_bin, Npoint, xc,yc,zc, r500, thickness)

Ncell = floor(Int32, (r500/(20*kpc)))
Ncell_box= 5
Nhalf_box=floor(Int32,Ncell_box/2)

box_deltad_2D, box_deltaT_2D, box_deltav_2D = box_filtering2D(d_2D, spec_T_2D, vx_2D, vy_2D, vz_2D, Ncell+Nhalf_box, Nhalf_box)


radial_deltad_2D = radial_deltad_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
radial_deltaT_2D = radial_deltaT_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
radial_deltav_2D = radial_deltav_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]

box_deltad_2D = box_deltad_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
box_deltaT_2D = box_deltaT_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
box_deltav_2D = box_deltav_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]

d_2D = d_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
spec_T_2D = spec_T_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]
v_2D = v_2D[xc-Ncell:xc+Ncell,yc-Ncell:yc+Ncell]

xc = floor(Int32, length(d_2D[:,1])/2 +1)
yc = floor(Int32, length(d_2D[1,:])/2 +1)


# Create the figure and subplots
fig, axes = subplots(3, 3, figsize=(12, 14), gridspec_kw=Dict("hspace" => 0.3, "wspace" => 0.0))

# First subplot
im1 = axes[1].pcolormesh(-Ncell*dx:dx:Ncell*dx, -Ncell*dx:dx:Ncell*dx, log10.(d_2D))
axes[1].axis("scaled")
axes[1].set_xlabel("y [Mpc]", fontsize=12)
axes[1].set_ylabel("z [Mpc]", fontsize=12)
axes[1].tick_params(axis="x", labelsize=12)
axes[1].tick_params(axis="y", labelsize=12)
divider = axgrid.make_axes_locatable(axes[1])
cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar1 = colorbar(im1, ax=axes[1], orientation="horizontal", cax=cax) # Adjust height as needed
cbar1.ax.xaxis.set_ticks_position("top")
cbar1.ax.xaxis.set_label_position("top")
cbar1.ax.tick_params(labelsize=10)
cbar1.set_label(L"Log (\rho \ [g/cm^3])", fontsize=12)
axes[1].set_xticklabels([])
axes[1].set_xlabel("")
axes[1].set_yticks([-0.5, 0.0, 0.5])
axes[1].set_yticklabels(["-0.5", "0.0", "0.5"])

# Second subplot
im2 = axes[2].pcolormesh(-Ncell*dx:dx:Ncell*dx, -Ncell*dx:dx:Ncell*dx, spec_T_2D*8.62e-8)
axes[2].axis("scaled")
axes[2].set_xlabel("y [Mpc]", fontsize=12)
axes[2].set_ylabel("z [Mpc]", fontsize=12)
divider = axgrid.make_axes_locatable(axes[2])
cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar2 = colorbar(im2, ax=axes[2], orientation="horizontal", cax=cax) # Adjust height as needed
cbar2.ax.xaxis.set_ticks_position("top")
cbar2.ax.xaxis.set_label_position("top")
cbar2.ax.tick_params(labelsize=10)
#cbar2.ax.yaxis.get_offset_text().set(size=13)
cbar2.set_label(L"T [keV]", fontsize=12)
axes[2].tick_params(axis="x", labelsize=12)
axes[2].tick_params(axis="y", labelsize=12)
axes[2].set_xticklabels([])
axes[2].set_xlabel("")
axes[2].set_yticks([-0.5, 0.0, 0.5])
axes[2].set_yticklabels(["-0.5", "0.0", "0.5"])

# Third subplot
im3 = axes[3].pcolormesh(-Ncell*dx:dx:Ncell*dx, -Ncell*dx:dx:Ncell*dx, log10.(v_2D*1.e-5))
axes[3].axis("scaled")
axes[3].set_xlabel("y [Mpc]", fontsize=12)
axes[3].set_ylabel("z [Mpc]", fontsize=12)
divider = axgrid.make_axes_locatable(axes[3])
cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar3 = colorbar(im3, ax=axes[3], orientation="horizontal", cax=cax) # Adjust height as needed
cbar3.ax.xaxis.set_ticks_position("top")
cbar3.ax.xaxis.set_label_position("top")
cbar3.ax.tick_params(labelsize=10)
cbar3.set_label(L"Log(v [km/s])", fontsize=12)
axes[3].tick_params(axis="x", labelsize=12)
axes[3].tick_params(axis="y", labelsize=12)
axes[3].set_yticks([-0.5, 0.0, 0.5])
axes[3].set_yticklabels(["-0.5", "0.0", "0.5"])

# Fourth subplot
im4 = axes[4].pcolormesh(-Ncell*dx:dx:Ncell*dx, -Ncell*dx:dx:Ncell*dx, radial_deltad_2D)
axes[4].axis("scaled")
axes[4].set_xlabel("y [Mpc]", fontsize=12)
axes[4].set_ylabel("z [Mpc]", fontsize=12)
divider = axgrid.make_axes_locatable(axes[4])
cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar4 = colorbar(im4, ax=axes[4], orientation="horizontal", cax=cax) # Adjust height as needed
cbar4.ax.xaxis.set_ticks_position("top")
cbar4.ax.xaxis.set_label_position("top")
cbar4.ax.tick_params(labelsize=10)
cbar4.set_label(L"\delta \rho/<\rho>", fontsize=12)
axes[4].tick_params(axis="x", labelsize=12)
axes[4].tick_params(axis="y", labelsize=12)
axes[4].set_yticklabels([])
axes[4].set_ylabel("")
axes[4].set_xticklabels([])
axes[4].set_xlabel("")


# Fifth subplot

im5 = axes[5].pcolormesh(-Ncell*dx:dx:Ncell*dx, -Ncell*dx:dx:Ncell*dx, radial_deltaT_2D)
axes[5].axis("scaled")
axes[5].set_xlabel("y [Mpc]", fontsize=12)
axes[5].set_ylabel("z [Mpc]", fontsize=12)
divider = axgrid.make_axes_locatable(axes[5])
cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar5 = colorbar(im5, ax=axes[5], orientation="horizontal", cax=cax) # Adjust height as needed
cbar5.ax.xaxis.set_ticks_position("top")
cbar5.ax.xaxis.set_label_position("top")
cbar5.ax.tick_params(labelsize=10)
cbar5.set_label(L"\delta T/<T>", fontsize=12)
axes[5].tick_params(axis="x", labelsize=12)
axes[5].tick_params(axis="y", labelsize=12)
axes[5].set_xticklabels([])
axes[5].set_xlabel("")
axes[5].set_yticklabels([])
axes[5].set_ylabel("")


# Sixth subplot

im6 = axes[6].pcolormesh(-Ncell*dx:dx:Ncell*dx, -Ncell*dx:dx:Ncell*dx, radial_deltav_2D)
axes[6].axis("scaled")
axes[6].set_xlabel("y [Mpc]", fontsize=12)
axes[6].set_ylabel("z [Mpc]", fontsize=12)
divider = axgrid.make_axes_locatable(axes[6])
cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar6 = colorbar(im6, ax=axes[6], orientation="horizontal", cax=cax) # Adjust height as needed
cbar6.ax.xaxis.set_ticks_position("top")
cbar6.ax.xaxis.set_label_position("top")
cbar6.ax.tick_params(labelsize=10)
cbar6.set_label(L"\delta v/<v>", fontsize=12)
axes[6].tick_params(axis="x", labelsize=12)
axes[6].tick_params(axis="y", labelsize=12)
axes[6].set_yticklabels([])
axes[6].set_ylabel("")


# Seventh subplot

im7 = axes[7].pcolormesh(-Ncell*dx:dx:Ncell*dx, -Ncell*dx:dx:Ncell*dx, box_deltad_2D)
axes[7].axis("scaled")
axes[7].set_xlabel("y [Mpc]", fontsize=12)
axes[7].set_ylabel("z [Mpc]", fontsize=12)
divider = axgrid.make_axes_locatable(axes[7])
cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar7 = colorbar(im7, ax=axes[7], orientation="horizontal", cax=cax) # Adjust height as needed
cbar7.ax.xaxis.set_ticks_position("top")
cbar7.ax.xaxis.set_label_position("top")
cbar7.ax.tick_params(labelsize=10)
cbar7.set_label(L"\delta \rho/<\rho>", fontsize=12)
axes[7].tick_params(axis="x", labelsize=12)
axes[7].tick_params(axis="y", labelsize=12)
axes[7].set_yticklabels([])
axes[7].set_ylabel("")
axes[7].set_xticklabels([])
axes[7].set_xlabel("")


# Eigth subplot

im8 = axes[8].pcolormesh(-Ncell*dx:dx:Ncell*dx, -Ncell*dx:dx:Ncell*dx, box_deltaT_2D)
axes[8].axis("scaled")
axes[8].set_xlabel("y [Mpc]", fontsize=12)
axes[8].set_ylabel("z [Mpc]", fontsize=12)
divider = axgrid.make_axes_locatable(axes[8])
cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar8 = colorbar(im8, ax=axes[8], orientation="horizontal", cax=cax) # Adjust height as needed
cbar8.ax.xaxis.set_ticks_position("top")
cbar8.ax.xaxis.set_label_position("top")
cbar8.ax.tick_params(labelsize=10)
cbar8.set_label(L"\delta T/<T>", fontsize=12)
axes[8].tick_params(axis="x", labelsize=12)
axes[8].tick_params(axis="y", labelsize=12)
axes[8].set_xticklabels([])
axes[8].set_xlabel("")
axes[8].set_yticklabels([])
axes[8].set_ylabel("")

# Nineth subplot

im9 = axes[9].pcolormesh(-Ncell*dx:dx:Ncell*dx, -Ncell*dx:dx:Ncell*dx, box_deltav_2D)
axes[9].axis("scaled")
axes[9].set_xlabel("y [Mpc]", fontsize=12)
axes[9].set_ylabel("z [Mpc]", fontsize=12)
divider = axgrid.make_axes_locatable(axes[9])
cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar9 = colorbar(im9, ax=axes[9], orientation="horizontal", cax=cax) # Adjust height as needed
cbar9.ax.xaxis.set_ticks_position("top")
cbar9.ax.xaxis.set_label_position("top")
cbar9.ax.tick_params(labelsize=10)
cbar9.set_label(L"\delta v/<v>", fontsize=12)
axes[9].tick_params(axis="x", labelsize=12)
axes[9].tick_params(axis="y", labelsize=12)
axes[9].set_yticklabels([])
axes[9].set_ylabel("")

show()
