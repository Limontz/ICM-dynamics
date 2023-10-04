using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using DataFrames
using LsqFit
using HypothesisTests
@pyimport mpl_toolkits.axes_grid1 as axgrid


main = "/home/marco/Scrivania/Ettori_project/"

# los = x
file1=string(main,"IT90_0/01_radial_filtering/new_x_IT90_0_2D_radial_filtering_fluctuations_radial_trend.txt")
file2=string(main,"IT90_1/01_radial_filtering/new_x_IT90_1_2D_radial_filtering_fluctuations_radial_trend.txt")
file3=string(main,"IT90_2/01_radial_filtering/new_x_IT90_2_2D_radial_filtering_fluctuations_radial_trend.txt")
file4=string(main,"IT90_3/01_radial_filtering/new_x_IT90_3_2D_radial_filtering_fluctuations_radial_trend.txt")
file5=string(main,"IT90_4/01_radial_filtering/new_x_IT90_4_2D_radial_filtering_fluctuations_radial_trend.txt")
file6=string(main,"IT92_0/01_radial_filtering/new_x_IT92_0_2D_radial_filtering_fluctuations_radial_trend.txt")
file7=string(main,"IT92_1/01_radial_filtering/new_x_IT92_1_2D_radial_filtering_fluctuations_radial_trend.txt")
file8=string(main,"IT92_2/01_radial_filtering/new_x_IT92_2_2D_radial_filtering_fluctuations_radial_trend.txt")

data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
data3=readdlm(file3,comments=true)
data4=readdlm(file4,comments=true)
data5=readdlm(file5,comments=true)
data6=readdlm(file6,comments=true)
data7=readdlm(file7,comments=true)
data8=readdlm(file8,comments=true)


density_x = cat(data8[:,2],cat(data7[:,2],cat(data6[:,2],cat(data5[:,2],cat(data4[:,2],cat(data3[:,2],cat(data2[:,2],data1[:,2],dims=(1,1)),dims=(1,1)),dims=(1,1)),
          dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
density2_x = cat(data8[:,3],cat(data7[:,3],cat(data6[:,3],cat(data5[:,3],cat(data4[:,3],cat(data3[:,3],cat(data2[:,3],data1[:,3],dims=(1,1)),dims=(1,1)),dims=(1,1)),
                dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
temperature_x = cat(data8[:,4],cat(data7[:,4],cat(data6[:,4],cat(data5[:,4],cat(data4[:,4],cat(data3[:,4],cat(data2[:,4],data1[:,4],dims=(1,1)),dims=(1,1)),
              dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity_x = cat(data8[:,5],cat(data7[:,5],cat(data6[:,5],cat(data5[:,5],cat(data4[:,5],cat(data3[:,5],cat(data2[:,5],data1[:,5],dims=(1,1)),dims=(1,1)),
           dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity2_x = cat(data8[:,6],cat(data7[:,6],cat(data6[:,6],cat(data5[:,6],cat(data4[:,6],cat(data3[:,6],cat(data2[:,6],data1[:,6],dims=(1,1)),dims=(1,1)),
           dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))



radius_x = cat(data8[:,1],cat(data7[:,1],cat(data6[:,1],cat(data5[:,1],cat(data4[:,1],cat(data3[:,1],cat(data2[:,1],data1[:,1],dims=(1,1)),dims=(1,1)),
               dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))

density_perturbed_x = cat(data8[:,2],cat(data7[:,2],cat(data5[:,2],cat(data4[:,2],data1[:,2],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
density2_perturbed_x = cat(data8[:,3],cat(data7[:,3],cat(data5[:,3],cat(data4[:,3],data1[:,3],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
temperature_perturbed_x = cat(data8[:,4],cat(data7[:,4],cat(data5[:,4],cat(data4[:,4],data1[:,4],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity_perturbed_x = cat(data8[:,5],cat(data7[:,5],cat(data5[:,5],cat(data4[:,5],data1[:,5],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity_perturbed2_x = cat(data8[:,6],cat(data7[:,6],cat(data5[:,6],cat(data4[:,6],data1[:,6],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
radius_perturbed_x = cat(data8[:,1],cat(data7[:,1],cat(data5[:,1],cat(data4[:,1],data1[:,1],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))

density_relaxed_x = cat(data6[:,2],cat(data3[:,2],data2[:,2],dims=(1,1)),dims=(1,1))
density2_relaxed_x = cat(data6[:,3],cat(data3[:,3],data2[:,3],dims=(1,1)),dims=(1,1))
temperature_relaxed_x = cat(data6[:,4],cat(data3[:,4],data2[:,4],dims=(1,1)),dims=(1,1))
velocity_relaxed_x = cat(data6[:,5],cat(data3[:,5],data2[:,5],dims=(1,1)),dims=(1,1))
velocity_relaxed2_x = cat(data6[:,6],cat(data3[:,6],data2[:,6],dims=(1,1)),dims=(1,1))
radius_relaxed_x = cat(data6[:,1],cat(data3[:,1],data2[:,1],dims=(1,1)),dims=(1,1))

wx = [8.9, 4.0, 5.3, 17.1, 11.5, 2.2, 13.6, 7.1]
wx = cat( fill(7.1,10), cat(fill(13.6,10),cat(fill(2.2,10),cat(fill(11.5,10), cat(fill(17.2,10),cat(fill(5.3,10),cat(fill(4.0,10),fill(8.9,10),
         dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1)),dims=(1,1)), dims=(1,1)), dims=(1,1) )


# los = y

file1=string(main,"IT90_0/01_radial_filtering/new_y_IT90_0_2D_radial_filtering_fluctuations_radial_trend.txt")
file2=string(main,"IT90_1/01_radial_filtering/new_y_IT90_1_2D_radial_filtering_fluctuations_radial_trend.txt")
file3=string(main,"IT90_2/01_radial_filtering/new_y_IT90_2_2D_radial_filtering_fluctuations_radial_trend.txt")
file4=string(main,"IT90_3/01_radial_filtering/new_y_IT90_3_2D_radial_filtering_fluctuations_radial_trend.txt")
file5=string(main,"IT90_4/01_radial_filtering/new_y_IT90_4_2D_radial_filtering_fluctuations_radial_trend.txt")
file6=string(main,"IT92_0/01_radial_filtering/new_y_IT92_0_2D_radial_filtering_fluctuations_radial_trend.txt")
file7=string(main,"IT92_1/01_radial_filtering/new_y_IT92_1_2D_radial_filtering_fluctuations_radial_trend.txt")
file8=string(main,"IT92_2/01_radial_filtering/new_y_IT92_2_2D_radial_filtering_fluctuations_radial_trend.txt")

data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
data3=readdlm(file3,comments=true)
data4=readdlm(file4,comments=true)
data5=readdlm(file5,comments=true)
data6=readdlm(file6,comments=true)
data7=readdlm(file7,comments=true)
data8=readdlm(file8,comments=true)


density_y = cat(data8[:,2],cat(data7[:,2],cat(data6[:,2],cat(data5[:,2],cat(data4[:,2],cat(data3[:,2],cat(data2[:,2],data1[:,2],dims=(1,1)),dims=(1,1)),dims=(1,1)),
          dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
density2_y = cat(data8[:,3],cat(data7[:,3],cat(data6[:,3],cat(data5[:,3],cat(data4[:,3],cat(data3[:,3],cat(data2[:,3],data1[:,3],dims=(1,1)),dims=(1,1)),dims=(1,1)),
                dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
temperature_y = cat(data8[:,4],cat(data7[:,4],cat(data6[:,4],cat(data5[:,4],cat(data4[:,4],cat(data3[:,4],cat(data2[:,4],data1[:,4],dims=(1,1)),dims=(1,1)),
              dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity_y = cat(data8[:,5],cat(data7[:,5],cat(data6[:,5],cat(data5[:,5],cat(data4[:,5],cat(data3[:,5],cat(data2[:,5],data1[:,5],dims=(1,1)),dims=(1,1)),
           dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity2_y = cat(data8[:,6],cat(data7[:,6],cat(data6[:,6],cat(data5[:,6],cat(data4[:,6],cat(data3[:,6],cat(data2[:,6],data1[:,6],dims=(1,1)),dims=(1,1)),
           dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))



radius_y = cat(data8[:,1],cat(data7[:,1],cat(data6[:,1],cat(data5[:,1],cat(data4[:,1],cat(data3[:,1],cat(data2[:,1],data1[:,1],dims=(1,1)),dims=(1,1)),
               dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))

density_perturbed_y = cat(data8[:,2],cat(data7[:,2],cat(data5[:,2],cat(data4[:,2],data1[:,2],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
density2_perturbed_y = cat(data8[:,3],cat(data7[:,3],cat(data5[:,3],cat(data4[:,3],data1[:,3],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
temperature_perturbed_y = cat(data8[:,4],cat(data7[:,4],cat(data5[:,4],cat(data4[:,4],data1[:,4],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity_perturbed_y = cat(data8[:,5],cat(data7[:,5],cat(data5[:,5],cat(data4[:,5],data1[:,5],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity_perturbed2_y = cat(data8[:,6],cat(data7[:,6],cat(data5[:,6],cat(data4[:,6],data1[:,6],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
radius_perturbed_y = cat(data8[:,1],cat(data7[:,1],cat(data5[:,1],cat(data4[:,1],data1[:,1],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))

density_relaxed_y = cat(data6[:,2],cat(data3[:,2],data2[:,2],dims=(1,1)),dims=(1,1))
density2_relaxed_y = cat(data6[:,3],cat(data3[:,3],data2[:,3],dims=(1,1)),dims=(1,1))
temperature_relaxed_y = cat(data6[:,4],cat(data3[:,4],data2[:,4],dims=(1,1)),dims=(1,1))
velocity_relaxed_y = cat(data6[:,5],cat(data3[:,5],data2[:,5],dims=(1,1)),dims=(1,1))
velocity_relaxed2_y = cat(data6[:,6],cat(data3[:,6],data2[:,6],dims=(1,1)),dims=(1,1))
radius_relaxed_y = cat(data6[:,1],cat(data3[:,1],data2[:,1],dims=(1,1)),dims=(1,1))

wy = [6.7, 10.3, 4.6, 4.1, 19.4, 3.9, 8.2, 5.2]
wy = cat( fill(5.2,10), cat(fill(8.2,10),cat(fill(3.9,10),cat(fill(19.4,10), cat(fill(4.1,10),cat(fill(4.6,10),cat(fill(10.3,10),fill(6.7,10),
         dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1)),dims=(1,1)), dims=(1,1)), dims=(1,1) )


#los=z


file1=string(main,"IT90_0/01_radial_filtering/new_z_IT90_0_2D_radial_filtering_fluctuations_radial_trend.txt")
file2=string(main,"IT90_1/01_radial_filtering/new_z_IT90_1_2D_radial_filtering_fluctuations_radial_trend.txt")
file3=string(main,"IT90_2/01_radial_filtering/new_z_IT90_2_2D_radial_filtering_fluctuations_radial_trend.txt")
file4=string(main,"IT90_3/01_radial_filtering/new_z_IT90_3_2D_radial_filtering_fluctuations_radial_trend.txt")
file5=string(main,"IT90_4/01_radial_filtering/new_z_IT90_4_2D_radial_filtering_fluctuations_radial_trend.txt")
file6=string(main,"IT92_0/01_radial_filtering/new_z_IT92_0_2D_radial_filtering_fluctuations_radial_trend.txt")
file7=string(main,"IT92_1/01_radial_filtering/new_z_IT92_1_2D_radial_filtering_fluctuations_radial_trend.txt")
file8=string(main,"IT92_2/01_radial_filtering/new_z_IT92_2_2D_radial_filtering_fluctuations_radial_trend.txt")

data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
data3=readdlm(file3,comments=true)
data4=readdlm(file4,comments=true)
data5=readdlm(file5,comments=true)
data6=readdlm(file6,comments=true)
data7=readdlm(file7,comments=true)
data8=readdlm(file8,comments=true)


density_z = cat(data8[:,2],cat(data7[:,2],cat(data6[:,2],cat(data5[:,2],cat(data4[:,2],cat(data3[:,2],cat(data2[:,2],data1[:,2],dims=(1,1)),dims=(1,1)),dims=(1,1)),
          dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
density2_z = cat(data8[:,3],cat(data7[:,3],cat(data6[:,3],cat(data5[:,3],cat(data4[:,3],cat(data3[:,3],cat(data2[:,3],data1[:,3],dims=(1,1)),dims=(1,1)),dims=(1,1)),
                dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
temperature_z = cat(data8[:,4],cat(data7[:,4],cat(data6[:,4],cat(data5[:,4],cat(data4[:,4],cat(data3[:,4],cat(data2[:,4],data1[:,4],dims=(1,1)),dims=(1,1)),
              dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity_z = cat(data8[:,5],cat(data7[:,5],cat(data6[:,5],cat(data5[:,5],cat(data4[:,5],cat(data3[:,5],cat(data2[:,5],data1[:,5],dims=(1,1)),dims=(1,1)),
           dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity2_z = cat(data8[:,6],cat(data7[:,6],cat(data6[:,6],cat(data5[:,6],cat(data4[:,6],cat(data3[:,6],cat(data2[:,6],data1[:,6],dims=(1,1)),dims=(1,1)),
           dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))



radius_z = cat(data8[:,1],cat(data7[:,1],cat(data6[:,1],cat(data5[:,1],cat(data4[:,1],cat(data3[:,1],cat(data2[:,1],data1[:,1],dims=(1,1)),dims=(1,1)),
               dims=(1,1)),dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))

density_perturbed_z = cat(data8[:,2],cat(data7[:,2],cat(data5[:,2],cat(data4[:,2],data1[:,2],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
density2_perturbed_z = cat(data8[:,3],cat(data7[:,3],cat(data5[:,3],cat(data4[:,3],data1[:,3],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
temperature_perturbed_z = cat(data8[:,4],cat(data7[:,4],cat(data5[:,4],cat(data4[:,4],data1[:,4],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity_perturbed_z = cat(data8[:,5],cat(data7[:,5],cat(data5[:,5],cat(data4[:,5],data1[:,5],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
velocity_perturbed2_z = cat(data8[:,6],cat(data7[:,6],cat(data5[:,6],cat(data4[:,6],data1[:,6],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))
radius_perturbed_z = cat(data8[:,1],cat(data7[:,1],cat(data5[:,1],cat(data4[:,1],data1[:,1],dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1))

density_relaxed_z = cat(data6[:,2],cat(data3[:,2],data2[:,2],dims=(1,1)),dims=(1,1))
density2_relaxed_z = cat(data6[:,3],cat(data3[:,3],data2[:,3],dims=(1,1)),dims=(1,1))
temperature_relaxed_z = cat(data6[:,4],cat(data3[:,4],data2[:,4],dims=(1,1)),dims=(1,1))
velocity_relaxed_z = cat(data6[:,5],cat(data3[:,5],data2[:,5],dims=(1,1)),dims=(1,1))
velocity_relaxed2_z = cat(data6[:,6],cat(data3[:,6],data2[:,6],dims=(1,1)),dims=(1,1))
radius_relaxed_z = cat(data6[:,1],cat(data3[:,1],data2[:,1],dims=(1,1)),dims=(1,1))

wz = [9.7, 5.6, 1.0, 15.4, 7.2, 5.4, 12.0, 7.7]
wz = cat( fill(7.7,10), cat(fill(12.0,10),cat(fill(5.4,10),cat(fill(7.2,10), cat(fill(15.4,10),cat(fill(1.0,10),cat(fill(5.6,10),fill(9.7,10),
         dims=(1,1)),dims=(1,1)),dims=(1,1)), dims=(1,1)),dims=(1,1)), dims=(1,1)), dims=(1,1) )



density = cat(density_x, cat(density_y, density_z, dims=(1,1)), dims=(1,1))
density2 = cat(density2_x, cat(density2_y, density2_z, dims=(1,1)), dims=(1,1))
temperature = cat(temperature_x, cat(temperature_y, temperature_z, dims=(1,1)), dims=(1,1))
velocity = cat(velocity_x, cat(velocity_y, velocity_z, dims=(1,1)), dims=(1,1))
velocity2 = cat(velocity2_x, cat(velocity2_y, velocity2_z, dims=(1,1)), dims=(1,1))
radius = cat(radius_x, cat(radius_y, radius_z, dims=(1,1)), dims=(1,1))

idx = findall(density2 .<= 2)
#density = density[idx]
density2 = density2[idx]
temperature_copy = temperature
temperature = temperature[idx]
#velocity = velocity[idx]
#velocity2 = velocity2[idx]


density_relaxed = cat(density_relaxed_x, cat(density_relaxed_y, density_relaxed_z, dims=(1,1)), dims=(1,1))
density2_relaxed = cat(density2_relaxed_x, cat(density2_relaxed_y, density2_relaxed_z, dims=(1,1)), dims=(1,1))
temperature_relaxed = cat(temperature_relaxed_x, cat(temperature_relaxed_y, temperature_relaxed_z, dims=(1,1)),dims=(1,1))
velocity_relaxed = cat(velocity_relaxed_x, cat(velocity_relaxed_y, velocity_relaxed_z, dims=(1,1)), dims=(1,1))
velocity_relaxed2 = cat(velocity_relaxed2_x, cat(velocity_relaxed2_y, velocity_relaxed2_z, dims=(1,1)), dims=(1,1))
radius_relaxed = cat(radius_relaxed_x, cat(radius_relaxed_y, radius_relaxed_z, dims=(1,1)), dims=(1,1))

density_perturbed = cat(density_perturbed_x, cat(density_perturbed_y, density_perturbed_z, dims=(1,1)), dims=(1,1))
density2_perturbed = cat(density2_perturbed_x, cat(density2_perturbed_y, density2_perturbed_z, dims=(1,1)), dims=(1,1))
temperature_perturbed = cat(temperature_perturbed_x, cat(temperature_perturbed_y, temperature_perturbed_z, dims=(1,1)), dims=(1,1))
velocity_perturbed = cat(velocity_perturbed_x, cat(velocity_perturbed_y, velocity_perturbed_z, dims=(1,1)), dims=(1,1))
velocity_perturbed2 = cat(velocity_perturbed2_x, cat(velocity_perturbed2_y, velocity_perturbed2_z, dims=(1,1)), dims=(1,1))
radius_perturbed = cat(radius_perturbed_x, cat(radius_perturbed_y, radius_perturbed_z, dims=(1,1)), dims=(1,1))

w = cat(wx, cat(wy, wz, dims=(1,1)), dims=(1,1))
w_copy = w
w_copy = w[idx]

m(x, p) =  p[1] .* x
p0=[1.]
fit = curve_fit(m, density2, temperature, p0)
param = fit.param[1]
#param = mean(temperature./density2)
sigma_par= stderror(fit)[1]
if sigma_par < 0.01
   sigma_par = 0.01
end
#res= fit.resid

R = Statistics.cor(density2, temperature)
p_value=pvalue(CorrelationTest(density2,temperature))
println("R= ", R)
println("p-value =", p_value)


println("----------- 2D FLUCTUATIONS: temperature-density -------------")
println("m = ", param, " +- ", sigma_par)

x = collect(minimum((density2)):0.01:maximum((density2)))
y = zeros(length(x))

for i in 1:length(x)
    y[i] = param*((x[i]))
end

y2 = zeros(length(x))

for i in 1:length(x)
    y2[i] = 0.19*((x[i]))
end

fig, ax1 = subplots(figsize=(6,7))
#fig, ax1 = subplots(figsize=(6,6))
fig.subplots_adjust(top=0.99,bottom=0.125,left=0.165, right=0.99)
#plot(density, temperature, linestyle="", marker=".", markersize="10.0",color="black")
scatter(density2, temperature, s=150, marker=".", c=w_copy, cmap="cividis")
#scatter(density_perturbed, temperature_perturbed, s=150, marker=".", facecolors="blue", edgecolors="black", label="Disturbed")
#scatter(density_relaxed, temperature_relaxed, s=150, marker=".", facecolors="orange", edgecolors="black", label="Relaxed")
ax=gca()
divider = axgrid.make_axes_locatable(ax)

plot(x,y,marker="", color="red")
ylabel(L"(\sigma_{T}/T(r))_{2D}", fontsize=20)
xlabel(L"(\sigma_{\rho^2}/\rho(r)^2)_{2D}", fontsize=20)
text(0.25, 0.37, string("m=0.20", L"\pm", round(sigma_par, digits=2)), fontsize= 18)
xticks([0.4, 0.8, 1.2, 1.6, 2.0], fontsize=18)
yticks([0.1, 0.2, 0.3, 0.4, 0.5], fontsize=18)
cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar=colorbar(orientation="horizontal", cax=cax)
cbar.ax.tick_params(labelsize=13)
cbar.ax.xaxis.set_ticks_position("top")
cbar.ax.xaxis.set_label_position("top")
##cbar.ax.yaxis.get_offset_text().set(size=20)
##cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"w[10^{-3}]",fontsize= 15)
##ylim(0,3)
#legend(fontsize=12)


show()
readline()
clf()



m(x, p) =  p[1] .* x
p0=[1.]
fit = curve_fit(m, density, velocity, p0)
param = fit.param[1]
#param = mean(temperature./density2)
sigma_par= stderror(fit)[1]
if sigma_par < 0.01
   sigma_par = 0.01
end
#res= fit.resid

#println("----------- 2D FLUCTUATIONS: temperature-density -------------")
println("m = ", param, " +- ", sigma_par)

x = collect(minimum((density)):0.01:1.0)
y = zeros(length(x))

for i in 1:length(x)
    y[i] = param*((x[i]))
end

R = Statistics.cor(density, velocity)
p_value=pvalue(CorrelationTest(density,velocity))
println("R= ", R)
println("p-value =", p_value)


fig, ax1 = subplots(figsize=(6,7))
#fig, ax1 = subplots(figsize=(6,6))
fig.subplots_adjust(top=0.99,bottom=0.125,left=0.165, right=0.99)
#plot(density, velocity, linestyle="", marker=".", markersize="10.0",color="black")
scatter(density, velocity, s=30, marker="s", c=w, cmap="cividis")
#scatter(density_perturbed, velocity_perturbed, s=30, marker="s", facecolors="blue", edgecolors="black", label="Disturbed")
#scatter(density_relaxed, velocity_relaxed, s=30, marker="s", facecolors="orange", edgecolors="black", label="Relaxed")

plot(x,y,marker="", color="red")
ylabel(L"(\sigma_{v}/c_s(r))_{2D}", fontsize=20)
xlabel(L"(\sigma_{\rho}/\rho(r))_{2D}", fontsize=20)
text(0.15, 0.47, string("m=0.50", L"\pm", round(sigma_par, digits=2)), fontsize= 18)
ax=gca()
divider = axgrid.make_axes_locatable(ax)
##ylim(0,3)
xticks([0.2, 0.4, 0.6, 0.8, 1.0], fontsize=18)
yticks([0.1, 0.2, 0.3, 0.4, 0.5], fontsize=18)
xlim(0.1,1.1)

cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar=colorbar(orientation="horizontal", cax=cax)
cbar.ax.tick_params(labelsize=13)
cbar.ax.xaxis.set_ticks_position("top")
cbar.ax.xaxis.set_label_position("top")
##cbar.ax.yaxis.get_offset_text().set(size=20)
##cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"w[10^{-3}]",fontsize= 15)

show()
readline()
clf()


m(x, p) =  p[1] .* x
p0=[1.]
fit = curve_fit(m, temperature_copy, velocity, p0)
param = fit.param[1]
#param = mean(temperature./density2)
sigma_par= stderror(fit)[1]
if sigma_par < 0.01
   sigma_par = 0.01
end
#res= fit.resid

println("----------- 2D FLUCTUATIONS: temperature-density -------------")
println("m = ", param, " +- ", sigma_par)

x = collect(minimum((temperature_copy)):0.01:maximum((temperature_copy)))
y = zeros(length(x))

for i in 1:length(x)
    y[i] = param[1]*((x[i]))
end


R = Statistics.cor(temperature_copy, velocity)
p_value=pvalue(CorrelationTest(temperature_copy, velocity))
println("R= ", R)
println("p-value =", p_value)


fig, ax1 = subplots(figsize=(6,7))
#fig, ax1 = subplots(figsize=(6,6))
fig.subplots_adjust(top=0.99,bottom=0.125,left=0.165, right=0.99)
#plot(temperature, velocity.^2, linestyle="", marker=".", markersize="10.0",color="black")
scatter(temperature_copy, velocity, s=30, marker="D", c=w, cmap="cividis")
#scatter(temperature_perturbed, velocity_perturbed2, s=30, marker="D", facecolors="blue", edgecolors="black", label="Disturbed")
#scatter(temperature_relaxed, velocity_relaxed2, s=30, marker="D", facecolors="orange", edgecolors="black", label="Relaxed")

plot(x,y,marker="", color="red")
ylabel(L"(\sigma_{v}/c_s(r)))_{2D}", fontsize=20)
xlabel(L"(\sigma_{T}/T(r))_{2D}", fontsize=20)
text(0.05, 0.45, string("m=", round(param, digits=2), L"\pm", round(sigma_par, digits=2)), fontsize= 18)

#ylim(0,3)
xticks([0.1, 0.2, 0.3, 0.4, 0.5], fontsize=18)
yticks([0.1, 0.2, 0.3, 0.4, 0.5], fontsize=18)
ax=gca()
divider = axgrid.make_axes_locatable(ax)
cax = divider[:append_axes]("top", size="5%", pad =0.05)
cbar=colorbar(orientation="horizontal", cax=cax)
cbar.ax.tick_params(labelsize=13)
cbar.ax.xaxis.set_ticks_position("top")
cbar.ax.xaxis.set_label_position("top")
##cbar.ax.yaxis.get_offset_text().set(size=20)
##cbar[:set_yticklabel](labelsize="large")
cbar[:set_label](L"w[10^{-3}]",fontsize= 15)

show()
