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


density_x = DataFrame( it900=data1[:,2], it901 = data2[:,2], it902 = data3[:,2], it903=data4[:,2], it904 = data5[:,2], it920 = data6[:,2], it921 = data7[:,2], it922=data8[:,2])
density2_x = DataFrame( it900=data1[:,3], it901 = data2[:,3], it902 = data3[:,3], it903=data4[:,3], it904 = data5[:,3], it920 = data6[:,3], it921 = data7[:,3], it922=data8[:,3])
temperature_x = DataFrame( it900=data1[:,4], it901 = data2[:,4], it902 = data3[:,4], it903=data4[:,4], it904 = data5[:,4], it920 = data6[:,4], it921 = data7[:,4], it922=data8[:,4])
velocity_x = DataFrame( it900=data1[:,5], it901 = data2[:,5], it902 = data3[:,5], it903=data4[:,5], it904 = data5[:,5], it920 = data6[:,5], it921 = data7[:,5], it922=data8[:,5])
velocity2_x = DataFrame( it900=data1[:,6], it901 = data2[:,6], it902 = data3[:,6], it903=data4[:,6], it904 = data5[:,6], it920 = data6[:,6], it921 = data7[:,6], it922=data8[:,6])




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


density_y = DataFrame( it900=data1[:,2], it901 = data2[:,2], it902 = data3[:,2], it903=data4[:,2], it904 = data5[:,2], it920 = data6[:,2], it921 = data7[:,2], it922=data8[:,2])
density2_y = DataFrame( it900=data1[:,3], it901 = data2[:,3], it902 = data3[:,3], it903=data4[:,3], it904 = data5[:,3], it920 = data6[:,3], it921 = data7[:,3], it922=data8[:,3])
temperature_y = DataFrame( it900=data1[:,4], it901 = data2[:,4], it902 = data3[:,4], it903=data4[:,4], it904 = data5[:,4], it920 = data6[:,4], it921 = data7[:,4], it922=data8[:,4])
velocity_y = DataFrame( it900=data1[:,5], it901 = data2[:,5], it902 = data3[:,5], it903=data4[:,5], it904 = data5[:,5], it920 = data6[:,5], it921 = data7[:,5], it922=data8[:,5])
velocity2_y = DataFrame( it900=data1[:,6], it901 = data2[:,6], it902 = data3[:,6], it903=data4[:,6], it904 = data5[:,6], it920 = data6[:,6], it921 = data7[:,6], it922=data8[:,6])






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


density_z = DataFrame( it900=data1[:,2], it901 = data2[:,2], it902 = data3[:,2], it903=data4[:,2], it904 = data5[:,2], it920 = data6[:,2], it921 = data7[:,2], it922=data8[:,2])
density2_z = DataFrame( it900=data1[:,3], it901 = data2[:,3], it902 = data3[:,3], it903=data4[:,3], it904 = data5[:,3], it920 = data6[:,3], it921 = data7[:,3], it922=data8[:,3])
temperature_z = DataFrame( it900=data1[:,4], it901 = data2[:,4], it902 = data3[:,4], it903=data4[:,4], it904 = data5[:,4], it920 = data6[:,4], it921 = data7[:,4], it922=data8[:,4])
velocity_z = DataFrame( it900=data1[:,5], it901 = data2[:,5], it902 = data3[:,5], it903=data4[:,5], it904 = data5[:,5], it920 = data6[:,5], it921 = data7[:,5], it922=data8[:,5])
velocity2_z = DataFrame( it900=data1[:,6], it901 = data2[:,6], it902 = data3[:,6], it903=data4[:,6], it904 = data5[:,6], it920 = data6[:,6], it921 = data7[:,6], it922=data8[:,6])


file1=string(main,"IT90_0/01_radial_filtering/new_IT90_0_3D_radial_filtering_fluctuations_radial_trend.txt")
file2=string(main,"IT90_1/01_radial_filtering/new_IT90_1_3D_radial_filtering_fluctuations_radial_trend.txt")
file3=string(main,"IT90_2/01_radial_filtering/new_IT90_2_3D_radial_filtering_fluctuations_radial_trend.txt")
file4=string(main,"IT90_3/01_radial_filtering/new_IT90_3_3D_radial_filtering_fluctuations_radial_trend.txt")
file5=string(main,"IT90_4/01_radial_filtering/new_IT90_4_3D_radial_filtering_fluctuations_radial_trend.txt")
file6=string(main,"IT92_0/01_radial_filtering/new_IT92_0_3D_radial_filtering_fluctuations_radial_trend.txt")
file7=string(main,"IT92_1/01_radial_filtering/new_IT92_1_3D_radial_filtering_fluctuations_radial_trend.txt")
file8=string(main,"IT92_2/01_radial_filtering/new_IT92_2_3D_radial_filtering_fluctuations_radial_trend.txt")

data1=readdlm(file1,comments=true)
data2=readdlm(file2,comments=true)
data3=readdlm(file3,comments=true)
data4=readdlm(file4,comments=true)
data5=readdlm(file5,comments=true)
data6=readdlm(file6,comments=true)
data7=readdlm(file7,comments=true)
data8=readdlm(file8,comments=true)


density = DataFrame( it900=data1[:,2], it901 = data2[:,2], it902 = data3[:,2], it903=data4[:,2], it904 = data5[:,2], it920 = data6[:,2], it921 = data7[:,2], it922=data8[:,2])
temperature = DataFrame( it900=data1[:,3], it901 = data2[:,3], it902 = data3[:,3], it903=data3[:,3], it904 = data5[:,3], it920 = data6[:,3], it921 = data7[:,3], it922=data8[:,3])
velocity = DataFrame( it900=data1[:,4], it901 = data2[:,4], it902 = data3[:,4], it903=data4[:,4], it904 = data4[:,4], it920 = data6[:,4], it921 = data7[:,4], it922=data8[:,4])
velocity2 = DataFrame( it900=data1[:,5], it901 = data2[:,5], it902 = data3[:,5], it903=data4[:,5], it904 = data5[:,5], it920 = data5[:,5], it921 = data7[:,5], it922=data8[:,5])

m(x, p) =  p[1] .* x
p0=[1.]

m_d_t_x = zeros(length(density_x[1,:]))
m_d_v_x = zeros(length(density_x[1,:]))
m_t_v_x = zeros(length(density_x[1,:]))
err_m_d_t_x = zeros(length(density_x[1,:]))
err_m_d_v_x = zeros(length(density_x[1,:]))
err_m_t_v_x = zeros(length(density_x[1,:]))

m_d_t_y = zeros(length(density_x[1,:]))
m_d_v_y = zeros(length(density_x[1,:]))
m_t_v_y = zeros(length(density_x[1,:]))
err_m_d_t_y = zeros(length(density_x[1,:]))
err_m_d_v_y = zeros(length(density_x[1,:]))
err_m_t_v_y = zeros(length(density_x[1,:]))

m_d_t_z = zeros(length(density_x[1,:]))
m_d_v_z = zeros(length(density_x[1,:]))
m_t_v_z = zeros(length(density_x[1,:]))
err_m_d_t_z = zeros(length(density_x[1,:]))
err_m_d_v_z = zeros(length(density_x[1,:]))
err_m_t_v_z = zeros(length(density_x[1,:]))

m_d_t = zeros(length(density_x[1,:]))
m_d_v = zeros(length(density_x[1,:]))
m_t_v = zeros(length(density_x[1,:]))
err_m_d_t = zeros(length(density_x[1,:]))
err_m_d_v = zeros(length(density_x[1,:]))
err_m_t_v = zeros(length(density_x[1,:]))

for i in 1:length(density_x[1,:])


    # x

    m_d_t_x[i] = mean(temperature_x[:,i]./density2_x[:,i])
    err_m_d_t_x[i] = std(temperature_x[:,i]./density2_x[:,i])

    m_d_v_x[i] = mean(velocity_x[:,i]./density2_x[:,i])
    err_m_d_v_x[i] = std(velocity_x[:,i]./density2_x[:,i])

    m_t_v_x[i] = mean(velocity_x[:,i]./temperature_x[:,i])
    err_m_t_v_x[i] = std(velocity_x[:,i]./temperature_x[:,i])

    # y
    m_d_t_y[i] = mean(temperature_y[:,i]./density2_y[:,i])
    err_m_d_t_y[i] = std(temperature_y[:,i]./density2_y[:,i])

    m_d_v_y[i] = mean(velocity_y[:,i]./density2_y[:,i])
    err_m_d_v_y[i] = std(velocity_y[:,i]./density2_y[:,i])

    m_t_v_y[i] = mean(velocity_y[:,i]./temperature_y[:,i])
    err_m_t_v_y[i] = std(velocity_y[:,i]./temperature_y[:,i])


    # z

    m_d_t_z[i] = mean(temperature_z[:,i]./density2_z[:,i])
    err_m_d_t_z[i] = std(temperature_z[:,i]./density2_z[:,i])

    m_d_v_z[i] = mean(velocity_z[:,i]./density2_z[:,i])
    err_m_d_v_z[i] = std(velocity_z[:,i]./density2_z[:,i])

    m_t_v_z[i] = mean(velocity_z[:,i]./temperature_z[:,i])
    err_m_t_v_z[i] = std(velocity_z[:,i]./temperature_z[:,i])

    #3D

    m_d_t[i] = mean(temperature[:,i]./density[:,i])
    err_m_d_t[i] = std(temperature[:,i]./density[:,i])

    m_d_v[i] = mean(velocity[:,i]./density[:,i])
    err_m_d_v[i] = std(velocity[:,i]./density[:,i])

    m_t_v[i] = mean(velocity[:,i]./temperature[:,i])
    err_m_t_v[i] = std(velocity[:,i]./temperature[:,i])


end



m_d_t_2D = zeros(length(m_d_v))
m_d_t_max = zeros(length(m_d_v))
m_d_t_min = zeros(length(m_d_v))

m_d_v_2D = zeros(length(m_d_v))
m_d_v_max = zeros(length(m_d_v))
m_d_v_min = zeros(length(m_d_v))

m_t_v_2D = zeros(length(m_d_v))
m_t_v_max = zeros(length(m_d_v))
m_t_v_min = zeros(length(m_d_v))


for i in 1:length(m_d_v)

    m_d_t_2D[i] = mean([m_d_t_x[i], m_d_t_y[i], m_d_t_z[i]])
    m_d_v_2D[i] = mean([m_d_v_x[i], m_d_v_y[i], m_d_v_z[i]])
    m_t_v_2D[i] = mean([m_t_v_x[i], m_t_v_y[i], m_d_t_z[i]])

    m_d_t_max[i] = maximum([m_d_t_x[i], m_d_t_y[i], m_d_t_z[i]])
    m_d_v_max[i] = maximum([m_d_v_x[i], m_d_v_y[i], m_d_v_z[i]])
    m_t_v_max[i] = maximum([m_t_v_x[i], m_t_v_y[i], m_d_t_z[i]])

    m_d_t_min[i] = minimum([m_d_t_x[i], m_d_t_y[i], m_d_t_z[i]])
    m_d_v_min[i] = minimum([m_d_v_x[i], m_d_v_y[i], m_d_v_z[i]])
    m_t_v_min[i] = minimum([m_t_v_x[i], m_t_v_y[i], m_d_t_z[i]])

end


l(x, p) =  p[1] .+ p[2].* x
p0=[0.5, 0]
fit = curve_fit(l, m_d_t, m_d_t_2D, p0)
q = fit.param[1]
slope = fit.param[2]
println(slope, " ", q)
x = collect(minimum((m_d_t)):0.01:maximum((m_d_t)))
y = zeros(length(x))
y2 = zeros(length(x))

for i in 1:length(x)
    y[i] = 0.5*((x[i])) + q
    y2[i] = slope*x[i] + q
end

index = range(1,8,step=1) |> collect
name = ["IT90_0","IT90_0", "IT90_1", "IT90_2", "IT90_3", "IT90_4", "IT92_0", "IT92_1", "IT92_2"]

fig, ax = subplots(figsize=(6,6))
fig.subplots_adjust(top=0.99,bottom=0.125,left=0.165, right=0.99)


errorbar(index.+0.1, m_d_t, yerr=err_m_d_t, marker=".", markersize=10.0, color="black", linestyle="", label = "3D")
errorbar(index, 2*m_d_t_2D, yerr=Array(2 .*[abs.(m_d_t_min - m_d_t_2D),2 .*abs.(m_d_t_max - m_d_t_2D)]), marker="s", markersize=6.0, color="green", linestyle="")
scatter(index, m_d_t_x, marker="D", s=50, facecolors="w", edgecolors="blue", label="2D")
scatter(index, m_d_t_y, marker="D", s=50, facecolors="w", edgecolors="blue")
scatter(index, m_d_t_z, marker="D", s=50, facecolors="w", edgecolors="blue")
xticks(fontsize=11)
yticks(fontsize=11)
ylabel(L"m_{2D}", fontsize = 11)
xlabel("Galaxy cluster", fontsize = 11)
ax.set_xticklabels(name)
plt.legend(fontsize=11)
show()
