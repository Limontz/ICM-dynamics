#profilo in k della potenza
function profile(d,temp, xc, yc, zc, r_max, rmin, Ncell_half, radial_bin)#(d_copy, temp, d_dm, vx, vy, vz, xc, yc, zc, G, kpc,
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
         #k=zeros(Nbox+1)

         P=zeros(Npoint, Npoint, Npoint)
         v_turb=zeros(Npoint, Npoint, Npoint)
         #S=zeros(Npoint, Npoint, Npoint)
         #v_r_mean = 0
         #v_t_mean = 0
         #M = 0


         @inbounds for i in xc - Ncell_half: xc + Ncell_half #(xc-N_half:xc+N_half) se voglio una boxcon Ncell centrata in xc,yc,zc)
             @inbounds for j in yc - Ncell_half: yc + Ncell_half
           @inbounds @simd for k in zc - Ncell_half: zc + Ncell_half

                     #@fastmath P[i ,j, k]=d_copy[i, j, k] * kb * temp[i, j, k] / (mu * mp) # compute pressure
                     #@fastmath S[i, j, k] = P[i,j,k] / d_copy[i,j,k]^(5. / 3.) # compute the modulus of the velocity
                     #@fastmath v_turb[i, j, k] = sqrt(vx[i, j, k]^2 + vy[i, j, k]^2 + vz[i, j, k]^2)
                     #@fastmath r=Int(floor(sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)))
                     @fastmath r=Int(floor(200*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)/r_max))
                               ra=convert(Int64,trunc(r))  # assign radial bin
                               #println(20*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)," ", 20*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)/r_max)

                               if ra <=1
                                  ra=1
                               end
                               if (ra >= radial_bin)
                                  ra=radial_bin
                               end

                               pk1[ra]+= d[i,j,k]
                               #pk2[ra]+= P[i,j,k]
                               #pk3[ra]+= S[i,j,k]
                               pk4[ra]+= 1.
                               pk5[ra]+= temp[i,j,k]


                               #pk7[ra]+= d_dm[i,j,k]
                               #pk8[ra]+= vx[i,j,k]^2
                              # v_t_mean += v_t[i,j,k]^2
                               #pk9[ra]+= beta[i,j,k]
                               #pk10[ra] += vy[i,j,k]^2
                               #v_r_mean += v_r[i,j,k]^2

                               #M += 1.



                               #if (d[i, j, k] >= min_d && d[i, j, k] <= max_d)
                                 # pk6[ra]+= v_turb[i,j,k]^2.0
                                #  N[ra]+=1.
                               #end
                           end
                       end
                   end


         for i in eachindex(pk1)
             pk1[i]=pk1[i]/pk4[i]
             #pk2[i]=pk2[i]/pk4[i]
             #pk3[i]=(pk3[i]/pk4[i])
             pk5[i]=pk5[i]/pk4[i]
             #pk6[i]= sqrt(pk6[i]/N[i])
             #pk7[i]= pk7[i]/pk4[i]
             #pk10[i]= (pk10[i]/pk4[i])
             #pk9[i]= pk9[i]/pk4[i]
         end


         #mass=zeros(radial_bin)
         #g=zeros(radial_bin)
         #dr=radius[3]-radius[2]
         #for i in eachindex(pk1)
        #     if(i==1)
                 #mass[i]=(4. *pi*(pk1[i]+pk7[i])*radius[i]^3.)/3. # total mass (dark matter + gas)
                 #g[i]=G*mass[i]/radius[i]^2.
             #else
            #     mass[i]=mass[i-1]+(4. *pi*(pk1[i]+pk7[i])*dr*radius[i]^2.) # total mass (daark matter + gas)
            #     g[i]=G*mass[i]/radius[i]^2. # "Newton's" gravitational acceleration
            # end
         #end


         return  pk1, pk5, radius
end
