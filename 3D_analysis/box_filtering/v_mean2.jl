#media della velocità pesata per la densità in box da 15 celle
function v_mean2(xc, yc, zc, Ncell, N_half, d, temp, vx, vy, vz, Npoint)

           num_x=0.
           num_y=0.
           num_z=0.
           den=0.

           num_d = 0.
           den_d = 0.
           num_t = 0.

           for i in 1:Ncell
               for j in 1:Ncell
                   for k in 1:Ncell

                       #print(i," ", j, " ", k)
                       num_x += d[xc-N_half+i,yc-N_half+j,zc-N_half+k]*vx[xc-N_half+i,yc-N_half+j,zc-N_half+k]
                       num_y += d[xc-N_half+i,yc-N_half+j,zc-N_half+k]*vy[xc-N_half+i,yc-N_half+j,zc-N_half+k]
                       num_z += d[xc-N_half+i,yc-N_half+j,zc-N_half+k]*vz[xc-N_half+i,yc-N_half+j,zc-N_half+k]
                       den += d[xc-N_half+i,yc-N_half+j,zc-N_half+k]
                      # println(d[xc-N_half+i-1,yc-N_half+i-1,zc-N_half+i-1])
                       num_d += d[xc-N_half+i,yc-N_half+j,zc-N_half+k]
                       num_t += temp[xc-N_half+i,yc-N_half+j,zc-N_half+k]
                       den_d += 1.

                   end
               end
           end
           vmean_x=num_x/den
           vmean_y=num_y/den
           vmean_z=num_z/den

           dmean = num_d / den_d
           tmean = num_t / den_d

           return(vmean_x,vmean_y,vmean_z, dmean, tmean)
end
