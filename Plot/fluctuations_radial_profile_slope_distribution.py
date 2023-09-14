import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 


from scipy.optimize import curve_fit

# func needs x and a as inputs
def func(x, a):
   return a*x

x = np.arange(0,0.5,0.02)
y = 0.33*x

popt, pcov = curve_fit(func, np.sqrt(x), y)

print(popt)
exit()



main = "/Users/marco/Desktop/The_Master/PhD/Ettori_project/"

file1=main + "Codes/Plot/2D_radial_filtering_fluctuation_radial_profile_slope_distribution.txt"
file2=main + "Codes/Plot/3D_radial_filtering_fluctuation_radial_profile_slope_distribution.txt"

name = ["IT90_0","IT90_0", "IT90_1", "IT90_2", "IT90_3", "IT90_4", "IT92_0", "IT92_1", "IT92_2"]
index = np.arange(1,9,1)
data=np.loadtxt(file1)
data3D=np.loadtxt(file2)

density_slope = data[:,0]
temperature_slope = data[:,1]
velocity_slope = data[:,2]

density_slope_3D = data3D[:,0]
temperature_slope_3D = data3D[:,1]
velocity_slope_3D = data3D[:,2]

perturbed_indices = [0,3,4,6,7]
relaxed_indices = [1,2,5]


d_mean_total_2D = np.mean(density_slope)
d_mean_perturbed_2D = np.mean(density_slope[perturbed_indices])
d_mean_relaxed_2D = np.mean(density_slope[relaxed_indices])

t_mean_total_2D = np.mean(temperature_slope)
t_mean_perturbed_2D = np.mean(temperature_slope[perturbed_indices])
t_mean_relaxed_2D = np.mean(temperature_slope[relaxed_indices])

v_mean_total_2D = np.mean(velocity_slope)
v_mean_perturbed_2D = np.mean(velocity_slope[perturbed_indices])
v_mean_relaxed_2D = np.mean(velocity_slope[relaxed_indices])



d_mean_total_3D = np.mean(density_slope_3D)
d_mean_perturbed_3D = np.mean(density_slope_3D[perturbed_indices])
d_mean_relaxed_3D = np.mean(density_slope_3D[relaxed_indices])

t_mean_total_3D = np.mean(temperature_slope_3D)
t_mean_perturbed_3D = np.mean(temperature_slope_3D[perturbed_indices])
t_mean_relaxed_3D = np.mean(temperature_slope_3D[relaxed_indices])

v_mean_total_3D = np.mean(velocity_slope_3D)
v_mean_perturbed_3D = np.mean(velocity_slope_3D[perturbed_indices])
v_mean_relaxed_3D = np.mean(velocity_slope_3D[relaxed_indices])



d_std_total_2D = np.std(density_slope)
d_std_perturbed_2D = np.std(density_slope[perturbed_indices])
d_std_relaxed_2D = np.std(density_slope[relaxed_indices])

t_std_total_2D = np.std(temperature_slope)
t_std_perturbed_2D = np.std(temperature_slope[perturbed_indices])
t_std_relaxed_2D = np.std(temperature_slope[relaxed_indices])

v_std_total_2D = np.std(velocity_slope)
v_std_perturbed_2D = np.std(velocity_slope[perturbed_indices])
v_std_relaxed_2D = np.std(velocity_slope[relaxed_indices])



d_std_total_3D = np.std(density_slope_3D)
d_std_perturbed_3D = np.std(density_slope_3D[perturbed_indices])
d_std_relaxed_3D = np.std(density_slope_3D[relaxed_indices])

t_std_total_3D = np.std(temperature_slope_3D)
t_std_perturbed_3D = np.std(temperature_slope_3D[perturbed_indices])
t_std_relaxed_3D = np.std(temperature_slope_3D[relaxed_indices])

v_std_total_3D = np.std(velocity_slope_3D)
v_std_perturbed_3D = np.std(velocity_slope_3D[perturbed_indices])
v_std_relaxed_3D = np.std(velocity_slope_3D[relaxed_indices])






print('2D slope:')
print('density slope (t,p,r)= ', round(d_mean_total_2D,2), '+-', round(d_std_total_2D,2), round(d_mean_perturbed_2D,2), '+-', round(d_std_perturbed_2D,2), round(d_mean_relaxed_2D,2), '+-', round(d_std_relaxed_2D,2))
print('temperature slope (t,p,r)= ', round(t_mean_total_2D,2), '+-', round(t_std_total_2D,2), round(t_mean_perturbed_2D,2), '+-', round(t_std_perturbed_2D,2), round(t_mean_relaxed_2D,2), '+-', round(d_std_relaxed_2D,2))
print('velocity slope (t,p,r)= ', round(v_mean_total_2D,2), '+-', round(v_std_total_2D,2),  round(v_mean_perturbed_2D,2), '+-', round(v_std_perturbed_2D,2), round(v_mean_relaxed_2D,2), '+-', round(v_std_relaxed_2D,2))

print('3D slope:')
print('density slope (t,p,r)= ', round(d_mean_total_3D,2), '+-', round(d_std_total_3D,2), round(d_mean_perturbed_3D,2), '+-', round(d_std_perturbed_3D,2), round(d_mean_relaxed_3D,2), '+-', round(d_std_relaxed_3D,2))
print('temperature slope (t,p,r)= ', round(t_mean_total_3D,2), '+-', round(t_std_total_3D,2), round(t_mean_perturbed_3D,2), '+-', round(t_std_perturbed_3D,2), round(t_mean_relaxed_3D,2), '+-', round(d_std_relaxed_3D,2))
print('velocity slope (t,p,r)= ', round(v_mean_total_3D,2), '+-', round(v_std_total_3D,2),  round(v_mean_perturbed_3D,2), '+-', round(v_std_perturbed_3D,2), round(v_mean_relaxed_3D,2), '+-', round(v_std_relaxed_3D,2))



fig = plt.figure(figsize = (6, 6))
ax = fig.add_subplot(2,1,1)
plt.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
plt.plot(index, density_slope, marker = '.', markersize = '10.0', linestyle='', color="black", label = "density")
plt.plot(index, temperature_slope, marker = 's', markersize = '6.0', linestyle='', color="blue", label = "temperature")
plt.plot(index, velocity_slope, marker = 'd', markersize = '6.0', linestyle='', color="red", label = "velocity")
plt.yticks(fontsize = 11)
plt.ylabel('slope(2D)', fontsize = 11)
plt.xticks(color='w')
plt.legend()

ax = fig.add_subplot(2,1,2)

plt.plot(index, density_slope_3D, marker = '.', markersize = '10.0', linestyle='', color="black", label = "density")
plt.plot(index, temperature_slope_3D, marker = 's', markersize = '6.0', linestyle='', color="blue", label = "temperature")
plt.plot(index, velocity_slope_3D, marker = 'd', markersize = '6.0', linestyle='', color="red", label = "velocity")
plt.xlabel('Galaxy cluster', fontsize=11)
plt.ylabel('slope(3D)', fontsize = 11)
plt.yticks(fontsize = 11)
ax.set_xticklabels(name)
plt.xticks(fontsize=11)


plt.show()





