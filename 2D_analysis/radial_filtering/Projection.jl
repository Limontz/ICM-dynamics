#Surface Brightness fluctuation

function projection(d, temp, vx, vy, vz, Npoint)

          #S_X = zeros(Npoint, Npoint, Npoint)
          d_2D = zeros(Npoint, Npoint)
          #S_X_2D = zeros(Npoint, Npoint)
          spec_T_2D = zeros(Npoint, Npoint)
          vx_2D = zeros(Npoint, Npoint)
          vy_2D = zeros(Npoint, Npoint)
          vz_2D = zeros(Npoint, Npoint)

@inbounds for i in 1:Npoint
@inbounds     for j in 1:Npoint
                  n = 0
                  den_T = 0
                  den_v = 0
@inbounds @simd   for k in 1:Npoint

                        #S_X[i, j, k] = (temp[i, j, k]^0.5) * d[i, j, k]^2

                        d_2D[i, j] += d[i, j, k]^2
                        n += 1
                        spec_T_2D[i,j] += d[i, j, k]^2 * temp[i,j,k]^(1/4)
                        den_T += d[i, j, k]^2 * temp[i,j,k]^(-3/4)
                        vx_2D[i,j] += d[i, j, k] * vx[i,j,k]
                        vy_2D[i,j] += d[i, j, k] * vy[i,j,k]
                        vz_2D[i,j] += d[i, j, k] * vz[i,j,k]
                        den_v += d[i, j, k]
                        #S_X_2D[j, k] += S_X[i, j, k]
                  end
                  d_2D[i, j] = sqrt(d_2D[i, j] / n)
                  spec_T_2D[i,j] = spec_T_2D[i,j] / den_T
                  vx_2D[i,j] = vx_2D[i,j] / den_v
                  vy_2D[i,j] = vy_2D[i,j] / den_v
                  vz_2D[i,j] = vz_2D[i,j] / den_v

              end

          end

          return(d_2D,spec_T_2D, vx_2D, vy_2D, vz_2D)
end
