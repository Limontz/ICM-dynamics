function filtering2D(d_2D, temp, vx, vy, vz, Npoint,  Nhalf_box)


    Ncell_box = floor(Int32,2*Nhalf_box)
    L = 2 * Nhalf_box * 20. * kpc
    lim = Nhalf_box

    T_turb_2D= zeros(length(d_2D[1,:]),length(d_2D[1,:]))
    d_turb_2D = zeros(length(d_2D[1,:]),length(d_2D[1,:]))
    v_turb_2D = zeros(length(d_2D[1,:]),length(d_2D[1,:]))
    lim=Nhalf_box

    @inbounds for i in lim:length(d_2D[1,:])-lim        #non cambia niente modificare il 26 perchè tanto le celle che uso
    @inbounds  @simd for j in lim:length(d_2D[1,:])-lim     # sono le celle centrali


                          kb = 1.38e-16
                          gamma = 5.0 / 3.0

                          d_mean, temp_mean, vx_mean, vy_mean, vz_mean  = v_mean2D(i, j, Ncell_box, Nhalf_box,d_2D,
                                                      temp, vx, vy, vz, Npoint) #evaluating the density/velocity mean on a box wide L

                          vx_turb = vx[i,j] - vx_mean
                          vy_turb = vy[i,j]-vy_mean
                          vz_turb = vz[i,j]-vz_mean

                          T_turb_2D[i,j] = (temp[i,j] - temp_mean) / temp_mean
                          d_turb_2D[i, j] = (d_2D[i, j] - d_mean) / d_mean

                          v_turb_2D[i,j] = sqrt( vx_turb^2. +  vy_turb^2. +  vz_turb^2.)/ sqrt((gamma * kb * temp[i,j] / (mu * mp)))
                          if (v_turb_2D[i, j] > 1.e3 || d_turb_2D[i, j] > 1.e3)
                              print(i," ",j)
                          end
                  end
              end

    return(d_turb_2D, T_turb_2D, v_turb_2D)

end
+
