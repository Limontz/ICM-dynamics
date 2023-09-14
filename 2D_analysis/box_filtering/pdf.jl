#Programma che fa genera una distribuzione
function pdf(dens, nr_bin)

    #nr_bin=floor(Int64,sqrt(length(d)))
    if (nr_bin ==1)
        nr_bin=3
    end
    bin=zeros(nr_bin)
    n=zeros(nr_bin)
    d_max, index = findmax(dens)
    d_min, index2 = findmin(dens)

    println(d_max,d_min)

@inbounds for i in 1:nr_bin
              bin[i]=d_min + (i-1.)*(d_max-d_min)/(nr_bin -1.) # generate a uniform grid of the denisty fluctuations
          end


          delta=bin[2] - bin[1]


          if (length(dens) > 1)
@inbounds @simd for i in eachindex(dens)
                    @fastmath ib=convert(Int64,trunc((dens[i]-d_min)/delta)) # find the bin in which the fluctuation falls
                    if(ib < 1)
                        ib=1
                    elseif(ib > nr_bin)
                        ib=nr_bin
                    end
                    n[ib]+= 1.
                end
            else
                n=1
            end

            nmax, ind=findmax(n)
            dmax=bin[ind]



    return(n,bin, dmax)
end
