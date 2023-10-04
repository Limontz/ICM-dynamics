#calcolo di v_turb per il profilo radiale del numero di Ri
function v_turb2(d, d_prof, vx_conv, vy_conv, vz_conv, temp, temp_prof, Npoint, radial_bin, xc, yc, zc)#, temp_prof)
                # vx_prof, vy_prof, vz_prof,temp_prof         #                v_prof
    mu = 0.62
    mp = 1.67e-24
    kb = 1.38e-16
    gamma = 5.0 / 3.0
    pc = 3.08e18
    kpc = pc * 1000
    #radial_bin=300
    Ncell= 1
    N_half=floor(Int32,Ncell/2)
    L= Ncell * 20. * kpc
    lim = 1
    vx_turb= zeros(Npoint,Npoint,Npoint)
    vy_turb= zeros(Npoint,Npoint,Npoint)
    vz_turb= zeros(Npoint,Npoint,Npoint)
    T_turb= zeros(Npoint,Npoint,Npoint)
    d_turb = zeros(Npoint, Npoint, Npoint)
    #d_turb_prof = zeros(Npoint, Npoint, Npoint)
    #log_d_turb = zeros(Npoint, Npoint, Npoint)
    #log_d_turb_prof = zeros(Npoint, Npoint, Npoint)
    v= zeros(Npoint,Npoint,Npoint)
   # v_box = zeros(Npoint-26,Npoint-26,Npoint-26)
@inbounds for i in lim:Npoint          #non cambia niente modificare il 26 perchè tanto le celle che uso
@inbounds     for j in lim:Npoint      # sono le celle centrali
@inbounds @simd   for k in lim:Npoint



                      kb = 1.38e-16
                      gamma = 5.0 / 3.0

            #@fastmath r=Int(floor(sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)))
                      @fastmath r=Int(floor(200*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)/r500))
                      ra=convert(Int64,trunc(r))  # assign radial bin
                      if ra <=1
                         ra=1
                      end
                      if ra >= radial_bin
                         ra=radial_bin
                      end

                      #vx_mean, vy_mean, vz_mean, d_mean, temp_mean = v_mean2(i, j, k, Ncell, N_half,d,
                                                                  #temp, vx_conv, vy_conv,
                                                                  #vz_conv, Npoint) #evaluating the density/velocity mean on a box wide L

                      #vx_turb[i,j,k] = (vx_conv[i,j,k] - vx_mean)# / temp_prof[ra]  #filtering each velocity component substracting the velocity mean to the original field

                      #vy_turb[i,j,k] = vy_conv[i,j,k]-vy_mean
                      #vz_turb[i,j,k] = vz_conv[i,j,k]-vz_mean

                      T_turb[i,j,k] = (temp[i,j,k] - temp_prof[ra])/temp_prof[ra]


                      #log_d_turb[i, j, k] = log(d[i, j, k]) - log(d_mean)
                      #log_d_turb_prof[k, j, i] = log(d[k, j, i]) - log(dprof[i])
                      #d_turb_prof[k, j, i] = (d[k, j, i] - dprof[i]) / (dprof[i])
                      d_turb[i, j, k] = (d[i, j, k] - d_prof[ra])/d_prof[ra]


                      vx_turb_prof = vx_conv[i,j,k] - vx_prof[ra] #filtering each velocity component substracting the velocity 3D profile
                      vy_turb_prof = vy_conv[i,j,k] - vy_prof[ra]
                      vz_turb_prof = vz_conv[i,j,k] - vz_prof[ra]

                      #v[i,j,k] = sqrt( vx_turb[i,j,k]^2. +  vy_turb[i,j,k]^2. +  vz_turb[i,j,k]^2.) / sqrt((gamma * kb * temp[i,j,k] / (mu * mp)))
                      vturb = sqrt( vx_turb_prof^2. + vy_turb_prof^2. + vz_turb_prof^2.) / sqrt((gamma * kb * temp[ra] / (mu * mp)))
                      #v_mean = sqrt( vx_mean^2. + vy_mean^2. +vz_mean^2.)
                      #v_box[i,j,k] = v[i, j, k] / sqrt((gamma * kb * temp_prof[ra] / (mu * mp)))  #evaluating the logarithmic velocity perturbations
                      #v_box_prof[i,j,k] = vturb[ra] / sqrt((gamma * kb * temp_prof[ra] / (mu * mp)))
                  end
              end
          end

    return(v, d_turb, T_turb, L)
end
