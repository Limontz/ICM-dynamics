function profile2D(d_2D, temp_2D, vx_2D, vy_2D, vz_2D, xc, yc, zc, r500, Ncell_half, radial_bin, thickness)

         pc = 3.08e18
         kpc = pc * 1000

         pk1 = zeros(radial_bin)
         pk2 = zeros(radial_bin)
         pk3 = zeros(radial_bin)
         pk4 = zeros(radial_bin)
         pk5 = zeros(radial_bin)
         N = zeros(radial_bin)
         println("radial bin = ",radial_bin)
         radius=zeros(radial_bin)
         #for i in 1:radial_bin
         #     radius[i]=rmin+(i-1.)*r_max/(radial_bin) # uniform grid. Each bin in connected to a radius value
         #end

         m = 0

@inbounds for i in xc - Ncell_half: xc + Ncell_half #(xc-N_half:xc+N_half) se voglio una boxcon Ncell centrata in xc,yc,zc)
@inbounds     for j in yc - Ncell_half: yc + Ncell_half

#@fastmath           r=Int(floor(sqrt((j - yc - 0.5)^2. +(i- xc -0.5)^2.)))
@fastmath           r=Int(floor(20*(1/thickness)*kpc*sqrt((j - yc - 0.5)^2. +(i- xc -0.5)^2.)/r500))
                    ra=convert(Int64,trunc(r))  # assign radial bin

                    if ra <=1
                       ra=1
                    end
                    if ra > radial_bin
                       #ra=radial_bin
                       @goto after
                    end

                    pk1[ra] += d_2D[i, j]
                    pk2[ra] += temp_2D[i, j]
                    pk3[ra] += vx_2D[i,j]
                    pk4[ra] += vy_2D[i,j]
                    pk5[ra] += vz_2D[i,j]
                    N[ra] += 1

                    @label after

                end
          end

          for i in eachindex(pk1)

             pk1[i] = pk1[i] / N[i]
             pk2[i] = pk2[i] / N[i]
             pk3[i] = pk3[i] / N[i]
             pk4[i] = pk4[i] / N[i]
             pk5[i] = pk5[i] / N[i]

          end



          return pk1, pk2, pk3, pk4, pk5, radius

end
