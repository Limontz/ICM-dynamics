#lettura e conversione
function read(snap, root)

         file1=string(root, "Ammassi/declust_z",snap)
         file2=string(root, "Ammassi/declust_v_z",snap)
         file3=string(root,"Ammassi/deDD",snap,".conv2")
         #file1 = string(root, "Ammassi/Fabrizio_cluster/test_profilo/", "dm_density.h5")
         #file2 = string(root, "Ammassi/Fabrizio_cluster/test_profilo/", "density.h5")
         #file3 = string(root, "Ammassi/Fabrizio_cluster/test_profilo/", "temperature.h5")

         conv=zeros(5)

         d_conv=h5read(file1, "Density")
         d_dm_conv=h5read(file1, "Dark_Matter_Density")
         temp=h5read(file1, "Temperature")
         #vx=h5read(file1,"vx")
         #vy=h5read(file1,"vy")
         #vz=h5read(file1,"vz")
         vx=h5read(file2,"x-velocity")
         vy=h5read(file2,"y-velocity")
         vz=h5read(file2,"z-velocity")

         #pressure = h5read(file1, "Pressure")

         conv=readdlm(file3)

         d_conv=d_conv*conv[4]
         d_dm_conv=d_dm_conv*conv[4]
         vx_conv=vx*conv[5]
         vy_conv=vy*conv[5]
         vz_conv=vz*conv[5]

         Npoint=length(d_conv[1,1,:])
         #…..Franco: in modo equivalente si può anche usare n=size(d)
         #………………………………………e Npoint=n[3] in questo caso

         return d_conv, d_dm_conv, temp, vx_conv, vy_conv, vz_conv, Npoint ,conv[4]#,d_dm_conv,temp,vx_conv,vy_conv,vz_conv,Npoint, conv[4]
end
