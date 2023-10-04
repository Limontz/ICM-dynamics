#calcolo di v_turb per il profilo radiale del numero di Ri
function v_turb2(d, vx_conv, vy_conv, vz_conv, temp, Npoint, Nhalf_box, xc, yc, zc)#, temp_prof)
                # vx_prof, vy_prof, vz_prof,temp_prof         #                v_prof
    mu = 0.62
    mp = 1.67e-24
    kb = 1.38e-16
    gamma = 5.0 / 3.0
    pc = 3.08e18
    kpc = pc * 1000


    Ncell_box = floor(Int32,2*Nhalf_box)
    L = 2 * Nhalf_box * 20. * kpc
    lim = Nhalf_box

    vx_turb= zeros(Npoint,Npoint,Npoint)
    vy_turb= zeros(Npoint,Npoint,Npoint)
    vz_turb= zeros(Npoint,Npoint,Npoint)
    T_turb= zeros(Npoint,Npoint,Npoint)
    d_turb = zeros(Npoint, Npoint, Npoint)
    v= zeros(Npoint,Npoint,Npoint)
   # v_box = zeros(Npoint-26,Npoint-26,Npoint-26)
@inbounds for i in lim:Npoint-lim          #non cambia niente modificare il 26 perchè tanto le celle che uso
@inbounds     for j in lim:Npoint-lim     # sono le celle centrali
@inbounds @simd   for k in lim:Npoint-lim

                      kb = 1.38e-16
                      gamma = 5.0 / 3.0


                      vx_mean, vy_mean, vz_mean, d_mean, temp_mean = v_mean2(i, j, k, Ncell_box, Nhalf_box,d,
                                                                  temp, vx_conv, vy_conv,
                                                                  vz_conv, Npoint) #evaluating the density/velocity mean on a box wide L

                      vx_turb[i,j,k] = (vx_conv[i,j,k] - vx_mean)# / temp_prof[ra]  #filtering each velocity component substracting the velocity mean to the original field
                      vy_turb[i,j,k] = vy_conv[i,j,k]-vy_mean
                      vz_turb[i,j,k] = vz_conv[i,j,k]-vz_mean

                      T_turb[i,j,k] = (temp[i,j,k] - temp_mean)/temp_mean

                      d_turb[i, j, k] = (d[i, j, k] - d_mean)/d_mean

                      v[i,j,k] = sqrt( vx_turb[i,j,k]^2. +  vy_turb[i,j,k]^2. +  vz_turb[i,j,k]^2.) / sqrt((gamma * kb * temp_mean / (mu * mp)))
                  end
              end
          end

    return(v, d_turb, T_turb, L)
end
