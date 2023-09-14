function mask2D(dens, temp, xc,yc, zc, lim, Npoint, radial_bin)

    new_dens = []
    new_temp = []

    for i in 1:Npoint
        for j in 1:Npoint

            r=Int(floor(200*kpc*sqrt((j - yc - 0.5)^2. +(i- xc -0.5)^2.)/r500))
            ra=convert(Int64,trunc(r))

            if ra <=1
               ra=1
            end
            if ra >= radial_bin
               ra=radial_bin
            end

            if (dens[i,j] < lim[ra])

                append!(new_dens, dens[i,j])
                append!(new_temp, temp[i,j])
            end
        end
    end

    return new_dens, new_temp

end
