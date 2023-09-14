function vmodule2D(vx,vy,vz,Npoint)

    v = zeros(Npoint,Npoint)

@inbounds for i in 1:Npoint
@inbounds @simd for j in 1:Npoint

                  v[i,j] = sqrt(vx[i,j]^2 + vy[i,j]^2 + vz[i,j]^2)

                end
          end

     return v

end
