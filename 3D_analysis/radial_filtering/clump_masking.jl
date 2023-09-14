#routine for clump masking

function masking(dens, Npoint, xc, yc, zc, radial_bin, r_max)


          #d_r, n_fluct = nfluct(Npoint, dens, radial_bin, xc, yc, zc, r_max)
          d_lim_95 = zeros(radial_bin)
          d_lim_99 = zeros(radial_bin)
          median = zeros(radial_bin)
          #nr_bin = 100

@inbounds for l in 1:radial_bin

              for_pdf = []

    @inbounds for k in 1:length(dens[1,1,:])
    @inbounds     for j in 1:length(dens[1,1,:])
   @inbounds @simd    for i in 1:length(dens[1,1,:])

                        #@fastmath v[i, j, k]=sqrt(vx_conv[i, j, k]^2.0 +vy_conv[i, j, k]^2.0 +vz_conv[i, j, k]^2.0)
                @fastmath r=Int(floor(200*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)/r_max))
                          ra=convert(Int64,trunc(r))
                          if(ra < 1)
                             ra=convert(Int64,1)
                           elseif(ra > radial_bin)
                                 ra=radial_bin
                           end

                           if(ra == l)

                              append!(for_pdf,dens[i,j,k])

                           end

                      end
                  end
              end

              n, bin, dmax = pdf(for_pdf, Int(floor(sqrt(length(for_pdf)))))

              median[l] = median!(for_pdf)
              m = 0
              n_tot = sum(n)

@inbounds     for i in 1:length(bin)

                   m += n[i]
                   if( m / n_tot >= 0.99)

                      d_lim_99[l] = (bin[i+1] + bin[i]) /2.
                      if (d_lim_99[l] < 5)
                         d_lim_99[l] = 5
                      end
                      break

                   end

              end
              m = 0
@inbounds     for i in 1:length(bin)

                     m += n[i]
                     if( m / n_tot >= 0.95)

                        d_lim_95[l] = (bin[i+1] + bin[i]) /2.
                        if (d_lim_95[l] < 5)
                           d_lim_95[l] = 5
                        end
                        break

                     end

              end

          end

          return(d_lim_95, d_lim_99, median)


end
