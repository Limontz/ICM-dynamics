function vmodule3D(vx,vy,vz,Npoint)

    v = zeros(Npoint,Npoint,Npoint)

@inbounds for i in 1:Npoint
@inbounds     for j in 1:Npoint
@inbounds @simd   for k in 1:Npoint

                      v[i,j,k] = sqrt(vx[i,j,k]^2 + vy[i,j,k]^2 + vz[i,j,k]^2)

                  end
              end
          end

    return v

end
