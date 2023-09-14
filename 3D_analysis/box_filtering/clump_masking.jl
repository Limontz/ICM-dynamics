#routine for clump masking

function masking(dens, Npoint, xc, yc, zc, radial_bin)

          d_r, n_fluct = nfluct(Npoint, dens, radial_bin, xc, yc, zc)
          d_lim_95 = zeros(radial_bin)
          d_lim_99 = zeros(radial_bin)
          median = zeros(radial_bin)
          #nr_bin = 100

@inbounds for k in 1:radial_bin

              n, bin, dmax = pdf(d_r, Int(floor(2*sqrt(length(d_r)))))
              println(bin)
              median[k] = median!(d_r)
              m = 0
              n_tot = sum(n)

@inbounds     for i in 1:length(bin)

                   m += n[i]
                   if( m / n_tot >= 0.99)

                      d_lim_99[k] = (bin[i+1] + bin[i]) /2.
                      break

                   end

              end
              m = 0
@inbounds     for i in 1:length(bin)

                     m += n[i]
                     if( m / n_tot >= 0.95)

                        d_lim_95[k] = bin[i+1]
                        break

                     end

              end

          end

          return(d_lim_95, d_lim_99, median)


end
