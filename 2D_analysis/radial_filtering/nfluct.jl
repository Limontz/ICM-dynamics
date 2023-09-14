#Calcolo del numero di fluttuazioni presenti ad un dato raggio

function nfluct(Npoint, d_conv, radial_bin, xc, yc, zc )

         d_r=zeros(radial_bin, convert(Int64,5e5))
         n_fluct=Array{Int64,1}(undef, radial_bin)
         #v=zeros(Npoint, Npoint, Npoint)

         for i in eachindex(n_fluct)
             n_fluct[i]=0
         end

         @inbounds for k in 1:Npoint
         @inbounds     for j in 1:Npoint
         @inbounds @simd    for i in 1:Npoint

                                #@fastmath v[i, j, k]=sqrt(vx_conv[i, j, k]^2.0 +vy_conv[i, j, k]^2.0 +vz_conv[i, j, k]^2.0)
                      @fastmath r=Int(floor(200*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)/r_max))
                                ra=convert(Int64,trunc(r))
                                if(ra < 1)
                                   ra=convert(Int64,1)
                                elseif(ra > radial_bin)
                                   ra=radial_bin
                                end


                                n_fluct[ra] += 1
                                n_fluct[ra]=floor(n_fluct[ra])

                                d_r[ra,n_fluct[ra]]=d[i,j,k]
                            end
                        end
                   end

                   return(d_r, n_fluct)
end
