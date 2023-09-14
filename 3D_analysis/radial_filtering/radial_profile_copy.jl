#profilo in k della potenza
function profile(d,temp, vx, vy, vz, xc, yc, zc, r_max, rmin, Ncell_half, radial_bin, thickness)#(d_copy, temp, d_dm, vx, vy, vz, xc, yc, zc, G, kpc,
                # r_max, rmin, Ncell, Ncell_half, max_d, min_d, radial_bin)

         G = 6.67e-8
         pc = 3.08e18
         kpc = pc * 1000
         kb = 1.38e-16
         mu = 0.62
         mp = 1.67e-24

         radius=zeros(radial_bin)
         for i in 1:radial_bin
             radius[i]=rmin+(i-1.)*r_max/(radial_bin) # uniform grid. Each bin in connected to a radius value
         end

         pk1=zeros(radial_bin)
         pk2=zeros(radial_bin)
         pk3=zeros(radial_bin)
         pk4=zeros(radial_bin)
         N=zeros(radial_bin)
         pk5=zeros(radial_bin)
         pk6=zeros(radial_bin)
         pk7= zeros(radial_bin)
         pk8= zeros(radial_bin)
         #pk9= zeros(radial_bin)
         pk10 = zeros(radial_bin)
         N = zeros(radial_bin)


         @inbounds for i in xc - Ncell_half: xc + Ncell_half #(xc-N_half:xc+N_half) se voglio una boxcon Ncell centrata in xc,yc,zc)
             @inbounds for j in yc - Ncell_half: yc + Ncell_half
           @inbounds       for k in zc - Ncell_half: zc + Ncell_half

                
                     @fastmath r=Int(floor(20*(1/thickness)*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)/r_max))
                               ra=convert(Int64,trunc(r))  # assign radial bin
                               #println(20*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)," ", 20*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)/r_max)

                               if ra <=1
                                  ra=1
                               end
                               if (ra >= radial_bin)
                                  #ra=radial_bin
                                  @goto skip
                               end

                               pk1[ra]+= d[i,j,k]
                               pk2[ra]+= vx[i,j,k]
                               pk3[ra]+= vy[i,j,k]
                               pk4[ra]+= 1.
                               pk5[ra]+= temp[i,j,k]
                               pk6[ra]+= vz[i,j,k]

                               @label skip
                           end
                       end
                   end


         for i in eachindex(pk1)
             pk1[i]=pk1[i]/pk4[i]
             pk2[i]=pk2[i]/pk4[i]
             pk3[i]=(pk3[i]/pk4[i])
             pk5[i]=pk5[i]/pk4[i]
             pk6[i]= (pk6[i]/pk4[i])
             #pk7[i]= pk7[i]/pk4[i]
             #pk10[i]= (pk10[i]/pk4[i])
             #pk9[i]= pk9[i]/pk4[i]
         end


         return  pk1, pk5, pk2, pk3, pk6, radius
end
