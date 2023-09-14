function v_mean2D(xc, yc, Ncell, N_half, d, temp, vx, vy, vz, Npoint)

           num_x=0.
           num_y=0.
           num_z=0.
           den=0.

           num_d = 0.
           den_d = 0.
           num_t = 0.

@inbounds  for i in 1:Ncell+1
@inbounds @simd for j in 1:Ncell+1

                       num_x += d[xc-N_half+i-1,yc-N_half+j-1]*vx[xc-N_half+i-1,yc-N_half+j-1]
                       num_y += d[xc-N_half+i-1,yc-N_half+j-1]*vy[xc-N_half+i-1,yc-N_half+j-1]
                       num_z += d[xc-N_half+i-1,yc-N_half+j-1]*vz[xc-N_half+i-1,yc-N_half+j-1]
                       den += d[xc-N_half+i-1,yc-N_half+j-1]
                      # println(d[xc-N_half+i-1,yc-N_half+i-1,zc-N_half+i-1])
                       num_d += d[xc-N_half+j-1,yc-N_half+j-1]
                       num_t += temp[xc-N_half+i-1,yc-N_half+j-1]
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
