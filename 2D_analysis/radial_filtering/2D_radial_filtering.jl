function filtering2D(d, d_prof, temp, temp_prof, vx, vx_prof, vy, vy_prof, vz, vz_prof, radial_bin, Npoint, xc,yc,zc, r500, thickness)

    T_turb_2D= zeros(Npoint,Npoint)
    d_turb_2D = zeros(Npoint,Npoint)
    v_turb_2D = zeros(Npoint,Npoint)
    lim=1
    @inbounds for i in lim:Npoint        #non cambia niente modificare il 26 perchè tanto le celle che uso
    @inbounds     for j in lim:Npoint     # sono le celle centrali


                          kb = 1.38e-16
                          gamma = 5.0 / 3.0

                          r=Int(floor(20*(1/thickness)*kpc*sqrt((j - yc - 0.5)^2. +(i- xc -0.5)^2.)/r500))
                          ra=convert(Int64,trunc(r))
                          if ra <=1
                             ra=1
                          end
                          if ra >= radial_bin
                             ra=radial_bin
                             #@goto after
                          end

                          vx_turb = vx[i,j]-vx_prof[ra]# / temp_prof[ra]  #filtering each velocity component substracting the velocity mean to the original field
                          vy_turb = vy[i,j]-vy_prof[ra]
                          vz_turb = vz[i,j]-vz_prof[ra]

                          T_turb_2D[i,j] = (temp[i,j] - temp_prof[ra]) / temp_prof[ra]
                          d_turb_2D[i, j] = (d[i, j] - d_prof[ra]) / d_prof[ra]
                          #println(i,j, " ", ra, " ", temp[i,j], " ", temp_prof[ra], " ", temp[i,j] - temp_prof[ra])

                          v_turb_2D[i,j] = sqrt( vx_turb^2. +  vy_turb^2. +  vz_turb^2.) / sqrt((gamma * kb * temp_prof[ra] / (mu * mp)))
                          #@label after
                  end
              end

    println("filtering fatto")

    return(d_turb_2D, T_turb_2D, v_turb_2D)

end
