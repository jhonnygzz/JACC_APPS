#function bhfp_threads_atomix(inputfile = get_input_filename_from_args();
#                             verbose = false)
#    (; ngauss, natom, xpnt, coef, geom) = parse_input_file(inputfile)
function bhfp_threads_atomix(data; verbose=false)
    (; ngauss, natom, xpnt, coef, geom) = data

    # this is needed because everytime coef and geom are modified in the function, it is modifying the original arrays from the benchmark_threads.jl!
    # (; ngauss, natom, xpnt, coef, geom) = data, creates references to the original arrays from benchmark_threads.jl
    coef = copy(coef) # make coef point to a new copy of coef, leaving the original untouched
    geom = copy(geom) # Nmake geom point to a new copy of geom, leaving the original untouched
        

    # println("dtol = ", BasicHFProxy.dtol)
    # #BasicHFProxy.dtol = 1.0e-12
    # println("rcut = ", BasicHFProxy.rcut)
    # println("tobohrs = ", BasicHFProxy.tobohrs)
    # println("sqrpi2 = ", BasicHFProxy.sqrpi2)   
    dens = Matrix{Float64}(undef, natom, natom)

    # fake ('random') density
    for i in 1:natom
        for j in 1:natom
            dens[i, j] = 0.1
        end
        dens[i, i] = 1.0
    end

    # normalize the primitive GTO weights
    for i in 1:ngauss
        coef[i] = coef[i] * (2.0 * xpnt[i])^0.75
    end

    # scale the geometry to Bohrs for energy calculations in AU
    for i in 1:natom
        geom[1, i] = geom[1, i] * tobohrs
        geom[2, i] = geom[2, i] * tobohrs
        geom[3, i] = geom[3, i] * tobohrs
    end

    fock = Matrix{Float64}(undef, natom, natom)
    for i in 1:natom
        for j in 1:natom
            fock[i, j] = 0.0
        end
    end

    # compute Schwarz Inequality factors for integral screening
    nn = ((natom^2) + natom) ÷ 2
    #nn = Int32((natom^2 + natom) ÷ 2)
    println("nn = ", nn)
    schwarz = Vector{Float64}(undef, nn)

    ij = 0
    for i in 1:natom
        for j in 1:i
            ij = ij + 1
            eri = ssss(i, j, i, j, ngauss, xpnt, coef, geom)
            schwarz[ij] = sqrt(abs(eri))
        end
    end

    t = _kernel_threaded_atomix!(fock, nn, schwarz, ngauss, xpnt, coef, geom, dens)

    # trace Fock with the density, print the 2e- energy
    erep = 0.0
    for i in 1:natom
        for j in 1:natom
            erep = erep + fock[i, j] * dens[i, j]
        end
    end
    E = erep * 0.5
    verbose && println("2e- energy= ", E)
    return E, t
end

function _kernel_threaded_atomix!(fock, nn, schwarz, ngauss, xpnt, coef, geom,
                                  dens)
    # The following loop (expanded to four indices, with permutational
    # symmetry) represents the kernel of Hartree-Fock calculations.
    # Integrals are screened to avoid small terms.
    #dtol = 1.0e-12
    nnnn = ((nn^2) + nn) ÷ 2
    #nnnn = Int32((nn^2 + nn) ÷ 2)
    println("nnnn= ", nnnn)

    start_time = time()
    Threads.@threads :static for ijkl in 1:nnnn
        # decompose triangular ijkl index into ij>=kl
        ij = isqrt(2 * ijkl)
        n = (ij * ij + ij) ÷ 2
        while (n < ijkl)
            ij = ij + 1
            n = (ij * ij + ij) ÷ 2
        end
        kl = ijkl - (ij * ij - ij) ÷ 2
        if schwarz[ij] * schwarz[kl] > dtol

            # decompose triangular ij index into i>=j
            i = isqrt(2 * ij)
            n = (i * i + i) ÷ 2
            while (n < ij)
                i = i + 1
                n = (i * i + i) ÷ 2
            end
            j = ij - (i * i - i) ÷ 2

            # decompose triangular kl index into k>=l
            k = isqrt(2 * kl)
            n = (k * k + k) ÷ 2
            while (n < kl)
                k = k + 1
                n = (k * k + k) ÷ 2
            end
            l = kl - (k * k - k) ÷ 2

            eri = 0.0

            for ib in 1:ngauss
                for jb in 1:ngauss
                    aij = 1.0 / (xpnt[ib] + xpnt[jb])
                    dij = coef[ib] * coef[jb] *
                          exp(-xpnt[ib] * xpnt[jb] * aij *
                              ((geom[1, i] - geom[1, j])^2 +
                               (geom[2, i] - geom[2, j])^2 +
                               (geom[3, i] - geom[3, j])^2)) * (aij^1.5)
                    if abs(dij) > dtol
                        xij = aij *
                              (xpnt[ib] * geom[1, i] + xpnt[jb] * geom[1, j])
                        yij = aij *
                              (xpnt[ib] * geom[2, i] + xpnt[jb] * geom[2, j])
                        zij = aij *
                              (xpnt[ib] * geom[3, i] + xpnt[jb] * geom[3, j])
                        for kb in 1:ngauss
                            for lb in 1:ngauss
                                akl = 1.0 / (xpnt[kb] + xpnt[lb])
                                dkl = dij * coef[kb] * coef[lb] *
                                      exp(-xpnt[kb] * xpnt[lb] * akl *
                                          ((geom[1, k] - geom[1, l])^2 +
                                           (geom[2, k] - geom[2, l])^2 +
                                           (geom[3, k] - geom[3, l])^2)) *
                                      (akl^1.5)
                                if (abs(dkl) > dtol)
                                    aijkl = (xpnt[ib] + xpnt[jb]) *
                                            (xpnt[kb] + xpnt[lb]) /
                                            (xpnt[ib] + xpnt[jb] + xpnt[kb] +
                                             xpnt[lb])
                                    tt = aijkl * ((xij -
                                           akl * (xpnt[kb] * geom[1, k] +
                                            xpnt[lb] * geom[1, l]))^2 +
                                          (yij -
                                           akl * (xpnt[kb] * geom[2, k] +
                                            xpnt[lb] * geom[2, l]))^2 +
                                          (zij -
                                           akl * (xpnt[kb] * geom[3, k] +
                                            xpnt[lb] * geom[3, l]))^2)
                                    f0t = sqrpi2
                                    if (tt > rcut)
                                        f0t = (tt^(-0.5)) *
                                              SpecialFunctions.erf(sqrt(tt))
                                    end
                                    eri = eri + dkl * f0t * sqrt(aijkl)
                                end
                            end
                        end
                    end
                end
            end

            if (i == j)
                eri = eri * 0.5
            end
            if (k == l)
                eri = eri * 0.5
            end
            if (i == k && j == l)
                eri = eri * 0.5
            end
            Atomix.@atomic fock[i, j] += dens[k, l] * eri * 4.0
            Atomix.@atomic fock[k, l] += dens[i, j] * eri * 4.0
            Atomix.@atomic fock[i, k] -= dens[j, l] * eri
            Atomix.@atomic fock[i, l] -= dens[j, k] * eri
            Atomix.@atomic fock[j, k] -= dens[i, l] * eri
            Atomix.@atomic fock[j, l] -= dens[i, k] * eri
            # @atomic fock[i, j] += dens[k, l] * eri * 4.0
            # @atomic fock[k, l] += dens[i, j] * eri * 4.0
            # @atomic fock[i, k] -= dens[j, l] * eri
            # @atomic fock[i, l] -= dens[j, k] * eri
            # @atomic fock[j, k] -= dens[i, l] * eri
            # @atomic fock[j, l] -= dens[i, k] * eri
        end
    end
    end_time = time()
    iteration_elapsed = end_time - start_time
    #println("iteration_elapsed ", iteration_elapsed)
    return iteration_elapsed
    #display(fock)
end