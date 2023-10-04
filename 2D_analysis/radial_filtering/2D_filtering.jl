function filtering2D(d, d_prof, temp, temp_prof, vx, vy, vz, radial_bin, Npoint, xc,yc,zc)

    Ncell= 5
    N_half=floor(Int32,Ncell/2)

    T_turb_2D= zeros(Npoint,Npoint)
    d_turb_2D = zeros(Npoint,Npoint)
    v_turb_2D = zeros(Npoint,Npoint)
    lim=1

    @inbounds for i in lim:Npoint        #non cambia niente modificare il 26 perchè tanto le celle che uso
    @inbounds  @simd for j in lim:Npoint     # sono le celle centrali


                          kb = 1.38e-16
                          gamma = 5.0 / 3.0

                @fastmath r=Int(floor(sqrt((j-yc -0.5)^2. +(i- xc -0.5)^2.)))
                          ra=convert(Int64,trunc(r))  # assign radial bin
                          if ra <=1
                             ra=1
                          end
                          if ra >= radial_bin
                             ra=radial_bin
                          end

                          d_mean, temp_mean, vx_mean, vy_mean, vz_mean  = v_mean2D(i, j, Ncell, N_half,d,
                                                      temp, vx, vy, vz, Npoint) #evaluating the density/velocity mean on a box wide L

                          vx_turb = vx[i,j] - vx_mean# / temp_prof[ra]  #filtering each velocity component substracting the velocity mean to the original field
                          vy_turb = vy[i,j]-vy_mean
                          vz_turb = vz[i,j]-vz_mean

                          T_turb_2D[i,j] = (temp[i,j] - temp_prof[ra]) / temp_prof[ra]
                          d_turb_2D[i, j] = (d[i, j] - d_prof[ra]) / d_prof[ra]
                          #println(i,j, " ", ra, " ", temp[i,j], " ", temp_prof[ra], " ", temp[i,j] - temp_prof[ra])



                          #vx_turb_prof = vx_conv[i,j,k] - vx_prof[ra] #filtering each velocity component substracting the velocity 3D profile
                          #vy_turb_prof = vy_conv[i,j,k] - vy_prof[ra]
                          #vz_turb_prof = vz_conv[i,j,k] - vz_prof[ra]
                          v = sqrt( vx_mean^2. +  vy_mean^2. +  vz_mean^2.)
                          cs = sqrt((gamma * kb * temp_prof[ra] / (mu * mp)))
                          v_turb_2D[i,j] = sqrt( vx_turb^2. +  vy_turb^2. +  vz_turb^2.) / cs

                  end
              end

    println("filtering fatto")

    return(d_turb_2D, T_turb_2D, v_turb_2D)

end
