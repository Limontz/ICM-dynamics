function v_mean2D(xc, yc, Ncell, N_half, d, temp, vx, vy, vz, Npoint)

           num_x=0.
           num_y=0.
           num_z=0.
           den=0.

           num_d = 0.
           den_d = 0.
           num_t = 0.

@inbounds  for i in 1:Ncell
@inbounds @simd for j in 1:Ncell

                       num_x += d[xc-N_half+i,yc-N_half+j]*vx[xc-N_half+i,yc-N_half+j]
                       num_y += d[xc-N_half+i,yc-N_half+j]*vy[xc-N_half+i,yc-N_half+j]
                       num_z += d[xc-N_half+i,yc-N_half+j]*vz[xc-N_half+i,yc-N_half+j]
                       den += d[xc-N_half+i,yc-N_half+j]
                      # println(d[xc-N_half+i-1,yc-N_half+i-1,zc-N_half+i-1])
                       num_d += d[xc-N_half+j,yc-N_half+j]
                       num_t += temp[xc-N_half+i,yc-N_half+j]
                       den_d += 1.

                   end
               end
           vmean_x=num_x/den
           vmean_y=num_y/den
           vmean_z=num_z/den

           dmean = num_d / den_d
           tmean = num_t / den_d

           return(dmean, tmean, vmean_x, vmean_y, vmean_z)
end
