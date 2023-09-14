function excision(d, temp, d_lim, median, xc, yc, zc, Npoint)

@inbounds for i in 1:Npoint
@inbounds     for j in 1:Npoint
@inbounds @simd   for k in 1:Npoint

            @fastmath r=Int(floor(200*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)/r_max))
                      ra=convert(Int64,trunc(r))
                      if(ra < 1)

                         ra=convert(Int64,1)

                      elseif(ra > radial_bin)

                         ra=radial_bin

                      end

                      if(d[i, j, k] > d_lim[ra])

                         d[i, j, k] = d_lim[ra]
                         #S_X[i, j, k] = (d[i, j, k]^2) * temp[i, j, k]^0.5

                      end

                   end
                end
             end

          return(d)

end
