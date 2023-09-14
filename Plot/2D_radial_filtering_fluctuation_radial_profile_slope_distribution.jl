using DelimitedFiles
using PyCall
using PyPlot
using Statistics
using DataFrames
using LsqFit

main = "/home/marco/Scrivania/Ettori_project/"

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

density = DataFrame(bin = data2[:,1], it900=data1[:,2],
                    it901 = data2[:,2],it902=data3[:,2],
                    it903=data4[:,2],it904=data5[:,2],
                    it920=data6[:,2],it921=data7[:,2],
                    it922=data8[:,2])

temperature = DataFrame(bin = data2[:,1], it900=data1[:,3],
                    it901 = data2[:,3],it902=data3[:,3],
                    it903=data4[:,3],it904=data5[:,3],
                    it920=data6[:,3],it921=data7[:,3],
                    it922=data8[:,3])

velocity = DataFrame(bin = data2[:,1], it900=data1[:,4],
                    it901 = data2[:,4],it902=data3[:,4],
                    it903=data4[:,4],it904=data5[:,4],
                    it920=data6[:,4],it921=data7[:,4],
                    it922=data8[:,4])


density_radial_slope = zeros(length(velocity[1,:])-1)
temperature_radial_slope = zeros(length(velocity[1,:])-1)
velocity_radial_slope = zeros(length(velocity[1,:])-1)

m(x, p) =  p[1] .+ p[2] .* x
p0=[0, 0.5]

for i in 1:length(density_radial_slope)

    fit = curve_fit(m, density[:,1]/10, density[:,i+1], p0)
    density_radial_slope[i] = fit.param[2]

    fit = curve_fit(m, temperature[:,1]/10, temperature[:,i+1], p0)
    temperature_radial_slope[i] = fit.param[2]

    fit = curve_fit(m, velocity[:,1]/10, velocity[:,i+1], p0)
    velocity_radial_slope[i] = fit.param[2]

end


output1=string("/home/marco/Scrivania/Ettori_project/Plots/2D_z_radial_filtering_fluctuation_radial_profile_slope_distribution.txt")
out1=open(output1,"w") do out1
writedlm( out1, [density_radial_slope temperature_radial_slope velocity_radial_slope])
end
close(output1)
