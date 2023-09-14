function total_clump_masking(d)

          nr_bin = 100

          n, bin, dmax = pdf(d, nr_bin)
          m = 0
          d_lim_95=0
          d_lim_99=0
          n_tot = sum(n)
@inbounds for i in 1:length(bin)

              m += n[i]
              if( m / n_tot >= 0.99)

                 d_lim_99 = (bin[i+1] + bin[i]) /2.
                 break

              end

          end
          m = 0
@inbounds for i in 1:length(bin)

               m += n[i]
               if( m / n_tot >= 0.80)

                  d_lim_95 = bin[i+1]
                  break

               end

          end

          return(n, bin, Int(floor(d_lim_95)), Int(floor(d_lim_99)))


end
