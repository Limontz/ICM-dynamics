function average_radial_perturb(dens, d_prof, temp, t_prof, vx, vx_prof, vy, vy_prof, vz, vz_prof, radial_bin, xc, yc, zc, rmax, thickness)

          mu = 0.62
          mp = 1.67e-24
          dens_shell_mean = zeros(radial_bin)
          temp_shell_mean = zeros(radial_bin)
          v_shell_mean = zeros(radial_bin)
          v_shell_mean2 = zeros(radial_bin)


@inbounds for l in 1:radial_bin

              dens_shell=[]
              temp_shell=[]
              v_shell=[]

    @inbounds for k in 1:length(dens[1,1,:])
    @inbounds     for j in 1:length(dens[1,1,:])
     @inbounds        for i in 1:length(dens[1,1,:])

                            #@fastmath v[i, j, k]=sqrt(vx_conv[i, j, k]^2.0 +vy_conv[i, j, k]^2.0 +vz_conv[i, j, k]^2.0)
                  @fastmath ra=Int(floor(20*(1/0.1)*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)/rmax))
                            ra=convert(Int64,trunc(ra))
                            if(ra < 1)
                              ra=convert(Int64,1)
                            elseif(ra > radial_bin)
                              #ra=radial_bin
                              @goto skip2
                            end

                            if(ra == l)

                               append!(dens_shell,dens[i,j,k])
                               append!(temp_shell,temp[i,j,k])
                               append!(v_shell,sqrt(vx[i,j,k]^2+vy[i,j,k]^2+vz[i,j,k]^2))

                            end

                            @label skip2

                        end
                    end
                end


                dens_shell_mean[l] = std(dens_shell)/d_prof[l]
                temp_shell_mean[l] = std(temp_shell)/t_prof[l]
                v_shell_mean[l] = std(v_shell)/ sqrt(gamma * kb * t_prof[l] / (mu * mp))
                v_shell_mean2[l] = std(v_shell.^2)/ (gamma * kb * t_prof[l] / (mu * mp))


          end

          return(dens_shell_mean, temp_shell_mean, v_shell_mean, v_shell_mean2)

end
