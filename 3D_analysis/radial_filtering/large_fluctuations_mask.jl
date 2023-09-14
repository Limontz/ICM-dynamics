function mask(dens, temp, v, xc,yc, zc, lim, radial_bin, r_max)

    new_dens = []
    new_temp = []
    new_v = []

    println(length(dens[:,1,1]))
    for i in 1:length(dens[:,1,1])
        for j in 1:length(dens[:,1,1])
            for k in 1:length(dens[:,1,1])

                r=Int(floor(200*kpc*sqrt((k -zc -0.5)^2. +(j-yc -0.5)^2. +(i- xc -0.5)^2.)/r_max))
                ra=convert(Int64,trunc(r))

                if ra <=1
                   ra=1
                end
                if (ra >= radial_bin)
                   ra=radial_bin
                end

                if (dens[i,j,k] < lim[ra])

                    append!(new_dens, dens[i,j,k])
                    append!(new_temp, temp[i,j,k])
                    append!(new_v, v[i,j,k])
                end

            end
        end
    end

    return new_dens, new_temp, new_v

end
