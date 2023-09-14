function average_radial_perturb(dens, temp, v, radial_bin, xc, yc, zc, rmax)

          dens_shell_mean = zeros(radial_bin)
          temp_shell_mean = zeros(radial_bin)
          v_shell_mean = zeros(radial_bin)
          param_dt = zeros(radial_bin)
          param_dv = zeros(radial_bin)
          param_tv = zeros(radial_bin)
          param_tv2 = zeros(radial_bin)
          #println(length(dens_shell_mean))

@inbounds for l in 1:radial_bin

              dens_shell=[]
              temp_shell=[]
              v_shell=[]

    @inbounds for i in 1:length(dens[1,:])
    @inbounds     for j in 1:length(dens[1,:])

                            #@fastmath v[i, j, k]=sqrt(vx_conv[i, j, k]^2.0 +vy_conv[i, j, k]^2.0 +vz_conv[i, j, k]^2.0)
                  @fastmath r=Int(floor(200*kpc*sqrt((j-yc -0.5)^2. +(i- xc -0.5)^2.)/rmax))
                            ra=convert(Int64,trunc(r))
                            if(ra < 1)
                              ra=convert(Int64,1)
                            elseif(ra > radial_bin)
                              #ra=radial_bin
                              @goto skip
                            end

                            if(ra == l)

                               append!(dens_shell,dens[i,j])
                               append!(temp_shell,temp[i,j])
                               append!(v_shell,v[i,j])

                            end
                            @label skip
                    end
                end

                m(x, p) = p[1] .* x

                p0=[1.]
                fit = curve_fit(m, ((dens_shell)), ((temp_shell)), p0)
                param_dt[l] = fit.param[1]

                p0=[1.]
                fit = curve_fit(m, (abs.(dens_shell)), (abs.(v_shell)), p0)
                param_dv[l]= fit.param[1]

                p0=[1.]
                fit = curve_fit(m, (abs.(temp_shell)), (abs.(v_shell)), p0)
                param_tv[l]= fit.param[1]

                p0=[1.]
                fit = curve_fit(m, (abs.(temp_shell)), (abs.(v_shell).^2), p0)
                param_tv2[l]= fit.param[1]

                dens_shell_mean[l] = std(dens_shell)
                temp_shell_mean[l] = std(temp_shell)
                v_shell_mean[l] = std(v_shell)

          end

          return(dens_shell_mean, temp_shell_mean, v_shell_mean, param_dt, param_dv, param_tv, param_tv2)

end
