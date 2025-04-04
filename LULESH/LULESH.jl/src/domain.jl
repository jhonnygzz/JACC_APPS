import JACC
import Pkg

# set backend
include("BasicLuleshPreferences.jl")
@static if endswith(BasicLuleshPreferences.backend, "cuda")
    # @TODO Julia Pkg.add will add target = :weakdeps in later versions
    # Pkg.add(; name = "CUDA", version = "v5.4.3")
    
    Pkg.add("CUDA")
    Pkg.update("CUDA")
    import CUDA
    # using CUDA
    # CUDA.set_runtime_version!(v"12.1.0")  # Replace with a supported CUDA version

    # pkg = Base.PkgId(Base.UUID("76a88914-d11a-5bdc-97e0-2f5a05c973a2"), "CUDA_Runtime_jll")
    # Base.compilecache(pkg)

    # CUDA.set_runtime_version!(v"11.2")  # Replace with your installed CUDA version
    println("Using CUDA as back end")
    println("CUDA version: ", CUDA.versioninfo())

    device = CUDA.device()
    println("device: $(CUDA.name(device))")

    devices = CUDA.devices()
    println("Available devices: ", devices)

    # CUDA.device!(1)

    current_device = CUDA.device()
    println("Current device: ", current_device)

elseif endswith(BasicLuleshPreferences.backend, "amdgpu")
    Pkg.add(; name = "AMDGPU", version = "v0.8.11")
    # Pkg.add(; name = "AMDGPU", version = "v1.0.3") 
    # Pkg.add("AMDGPU")
    # Pkg.update("AMDGPU")
    import AMDGPU
    # using AMDGPU
    println("Using AMDGPU as back end")
    # devices = AMDGPU.devices()

    
    # device = devices[2]
    # AMDGPU.device!(device)
    # # device = AMDGPU.device()
    # println("device: $device")

    # # Check available devices
    # devices = AMDGPU.devices()
    # println("Number of devices: ", length(devices))

    # if !isempty(devices)
    #     # Print details about the first device
    #     device = devices[1]
    #     println("Device 1:")
    #     println("  Device type: ", typeof(device))
    #     println("  Device name: ", device)
    # end

elseif endswith(BasicLuleshPreferences.backend, "threads")
    using Base.Threads
    println("Using threads as back end")
    # Threads.nthreads() = 192
    # println("Number of threads: ", Threads.nthreads())
end



# access elements for comms
get_delv_xi(idx::IndexT, dom::AbstractDomain) = dom.d_delv_xi[idx]
get_delv_eta(idx::IndexT, dom::AbstractDomain) = dom.d_delv_eta[idx]
get_delv_zeta(idx::IndexT, dom::AbstractDomain) = dom.d_delv_zeta[idx]

get_x(idx::IndexT, dom::AbstractDomain) = dom.d_x[idx]
get_y(idx::IndexT, dom::AbstractDomain) = dom.d_y[idx]
get_z(idx::IndexT, dom::AbstractDomain) = dom.d_z[idx]

get_xd(idx::IndexT, dom::AbstractDomain) = dom.d_xd[idx]
get_yd(idx::IndexT, dom::AbstractDomain) = dom.d_yd[idx]
get_zd(idx::IndexT, dom::AbstractDomain) = dom.d_zd[idx]

get_fx(idx::IndexT, dom::AbstractDomain) = dom.d_fx[idx]
get_fy(idx::IndexT, dom::AbstractDomain) = dom.d_fy[idx]
get_fz(idx::IndexT, dom::AbstractDomain) = dom.d_fz[idx]

# assume communication to 6 neighbors by default
m_rowMin(domain::Domain) = (domain.m_rowLoc == 0)             ? false : true
m_rowMax(domain::Domain) = (domain.m_rowLoc == domain.m_tp-1) ? false : true
m_colMin(domain::Domain) = (domain.m_colLoc == 0)             ? false : true
m_colMax(domain::Domain) = (domain.m_colLoc == domain.m_tp-1) ? false : true
m_planeMin(domain::Domain) = (domain.m_planeLoc == 0)         ? false : true
m_planeMax(domain::Domain) = (domain.m_planeLoc == domain.m_tp-1) ? false : true

# host access
get_nodalMass(idx::IndexT, dom::AbstractDomain) = dom.nodalMass[idx]

colLoc(dom::AbstractDomain) = dom.m_colLoc
rowLoc(dom::AbstractDomain) = dom.m_rowLoc
planeLoc(dom::AbstractDomain) = dom.m_planeLoc
tp(dom::AbstractDomain) = dom.m_tp

function allocateNodalPersistent!(domain, domNodes)
    fill!(resize!(domain.x, domNodes),0)   # coordinates
    fill!(resize!(domain.y, domNodes),0)
    fill!(resize!(domain.z, domNodes),0)

    fill!(resize!(domain.xd, domNodes),0)  # velocities
    fill!(resize!(domain.yd, domNodes),0)
    fill!(resize!(domain.zd, domNodes),0)

    fill!(resize!(domain.xdd, domNodes),0) # accelerations
    fill!(resize!(domain.ydd, domNodes),0) # accelerations
    fill!(resize!(domain.zdd, domNodes),0) # accelerations

    fill!(resize!(domain.fx, domNodes),0)   # forces
    fill!(resize!(domain.fy, domNodes),0)
    fill!(resize!(domain.fz, domNodes),0)

    fill!(resize!(domain.dfx, domNodes),0)  # AD derivative of the forces
    fill!(resize!(domain.dfy, domNodes),0)
    fill!(resize!(domain.dfz, domNodes),0)

    fill!(resize!(domain.nodalMass, domNodes),0)  # mass
    return nothing
end

function allocateElemPersistent!(domain, domElems)
    resize!(domain.matElemlist, domElems)  # material indexset
    resize!(domain.nodelist, 8*domElems)   # elemToNode connectivity
    fill!(domain.nodelist, 0)

    resize!(domain.lxim, domElems)  # elem connectivity through face g
    resize!(domain.lxip, domElems)
    resize!(domain.letam, domElems)
    resize!(domain.letap, domElems)
    resize!(domain.lzetam, domElems)
    resize!(domain.lzetap, domElems)

    resize!(domain.elemBC, domElems)   # elem face symm/free-surf flag g

    resize!(domain.e, domElems)    # energy g
    resize!(domain.p, domElems)    # pressure g

    resize!(domain.d_e, domElems)  # AD derivative of energy E g

    resize!(domain.q, domElems)    # q g
    resize!(domain.ql, domElems)   # linear term for q g
    resize!(domain.qq, domElems)   # quadratic term for q g
    resize!(domain.v, domElems)      # relative volume g

    resize!(domain.volo, domElems)   # reference volume g
    resize!(domain.delv, domElems)   # m_vnew - m_v g
    resize!(domain.vdov, domElems)   # volume derivative over volume g

    resize!(domain.arealg, domElems)   # elem characteristic length g

    resize!(domain.ss, domElems)       # "sound speed" g

    resize!(domain.elemMass, domElems)   # mass g
    return nothing
end

function initializeFields!(domain)
    # Basic Field Initialization

    fill!(domain.ss,0.0);
    fill!(domain.e,0.0)
    fill!(domain.p,0.0)
    fill!(domain.q,0.0)
    fill!(domain.v,1.0)

    fill!(domain.d_e,0.0)

    fill!(domain.xd,0.0)
    fill!(domain.yd,0.0)
    fill!(domain.zd,0.0)

    fill!(domain.xdd,0.0)
    fill!(domain.ydd,0.0)
    fill!(domain.zdd,0.0)

    fill!(domain.nodalMass,0.0)
end

function buildMesh!(domain, nx, edgeNodes, edgeElems, domNodes, domElems, x, y, z, nodelist)
    meshEdgeElems = domain.m_tp*nx

    resize!(x, domNodes)
    resize!(y, domNodes)
    resize!(z, domNodes)
    # initialize nodal coordinates
    # INDEXING
    nidx::IndexT = 1
    tz = 1.125*(domain.m_planeLoc*nx)/meshEdgeElems
    for plane in 1:edgeNodes
        ty = 1.125*(domain.m_rowLoc*nx)/meshEdgeElems
        for row in 1:edgeNodes
        tx = 1.125*(domain.m_colLoc*nx)/meshEdgeElems
            for col in 1:edgeNodes
                x[nidx] = tx
                y[nidx] = ty
                z[nidx] = tz
                nidx+=1
                # tx += ds ; // may accumulate roundoff...
                tx = 1.125*(domain.m_colLoc*nx+col)/meshEdgeElems
            end
        #// ty += ds ;  // may accumulate roundoff...
        ty = 1.125*(domain.m_rowLoc*nx+row)/meshEdgeElems
        end
        #// tz += ds ;  // may accumulate roundoff...
        tz = 1.125*(domain.m_planeLoc*nx+plane)/meshEdgeElems
    end

    copyto!(domain.x, x)
    copyto!(domain.y, y)
    copyto!(domain.z, z)
    resize!(nodelist, domElems*8);

    # embed hexehedral elements in nodal point lattice
    # INDEXING
    zidx::IndexT = 0
    nidx = 0
    for plane in 1:edgeElems
        for row in 1:edgeElems
            for col in 1:edgeElems
                nodelist[8*zidx+1] = nidx
                nodelist[8*zidx+2] = nidx                                   + 1
                nodelist[8*zidx+3] = nidx                       + edgeNodes + 1
                nodelist[8*zidx+4] = nidx                       + edgeNodes
                nodelist[8*zidx+5] = nidx + edgeNodes*edgeNodes
                nodelist[8*zidx+6] = nidx + edgeNodes*edgeNodes             + 1
                nodelist[8*zidx+7] = nidx + edgeNodes*edgeNodes + edgeNodes + 1
                nodelist[8*zidx+8] = nidx + edgeNodes*edgeNodes + edgeNodes
                zidx+=1
                nidx+=1
            end
        nidx+=1
        end
    nidx+=edgeNodes
    end
    copyto!(domain.nodelist, nodelist)
end

function setupConnectivityBC!(domain::Domain, edgeElems)
    domElems = domain.numElem;

    lxim = Vector{IndexT}(undef, domElems)
    lxip = Vector{IndexT}(undef, domElems)
    letam = Vector{IndexT}(undef, domElems)
    letap = Vector{IndexT}(undef, domElems)
    lzetam = Vector{IndexT}(undef, domElems)
    lzetap = Vector{IndexT}(undef, domElems)

    # set up elemement connectivity information
    lxim[1] = 0 ;
    for i in 2:domElems
       lxim[i]   = i-2
       lxip[i-1] = i-1
    end
    # MAYBE
    lxip[domElems] = domElems-1

    # INDEXING
    for i in 1:edgeElems
       letam[i] = i-1
       letap[domElems-edgeElems+i] = domElems-edgeElems+i-1
    end

    for i in (edgeElems+1):domElems
       letam[i] = i-edgeElems-1
       letap[i-edgeElems] = i-1
    end

    for i in 1:edgeElems*edgeElems
       lzetam[i] = i-1
       lzetap[domElems-edgeElems*edgeElems+i] = domElems-edgeElems*edgeElems+i-1
    end

    for i in (edgeElems*edgeElems+1):domElems
       lzetam[i] = i - edgeElems*edgeElems-1
       lzetap[i-edgeElems*edgeElems] = i-1
    end


    # set up boundary condition information
    elemBC = Vector{IndexT}(undef, domElems)
    for i in 1:domElems
        elemBC[i] = 0   # clear BCs by default
    end

    ghostIdx = [typemin(IndexT) for i in 1:6]::Vector{IndexT} # offsets to ghost locations

    pidx = domElems
    if m_planeMin(domain) != 0
        ghostIdx[1] = pidx
        pidx += domain.sizeX*domain.sizeY
    end

    if m_planeMax(domain) != 0
        ghostIdx[2] = pidx
        pidx += domain.sizeX*domain.sizeY
    end

    if m_rowMin(domain) != 0
        ghostIdx[3] = pidx
        pidx += domain.sizeX*domain.sizeZ
    end

    if m_rowMax(domain) != 0
        ghostIdx[4] = pidx
        pidx += domain.sizeX*domain.sizeZ
    end

    if m_colMin(domain) != 0
        ghostIdx[5] = pidx
        pidx += domain.sizeY*domain.sizeZ
    end

    if m_colMax(domain) != 0
        ghostIdx[6] = pidx
    end

    # symmetry plane or free surface BCs
    for i in 1:edgeElems
        planeInc = (i-1)*edgeElems*edgeElems
        rowInc   = (i-1)*edgeElems
        for j in 1:edgeElems
            if domain.m_planeLoc == 0
                elemBC[rowInc+j] |= ZETA_M_SYMM
            else
                elemBC[rowInc+j] |= ZETA_M_COMM
                lzetam[rowInc+j] = ghostIdx[1] + rowInc + (j-1)
            end

            if domain.m_planeLoc == domain.m_tp-1
                elemBC[rowInc+j+domElems-edgeElems*edgeElems] |= ZETA_P_FREE
            else
                elemBC[rowInc+j+domElems-edgeElems*edgeElems] |= ZETA_P_COMM
                lzetap[rowInc+j+domElems-edgeElems*edgeElems] = ghostIdx[2] + rowInc + (j-1)
            end

            if domain.m_rowLoc == 0
                elemBC[planeInc+j] |= ETA_M_SYMM
            else
                elemBC[planeInc+j] |= ETA_M_COMM
                letam[planeInc+j] = ghostIdx[3] + rowInc + (j-1)
            end

            if domain.m_rowLoc == domain.m_tp-1
                elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_FREE
            else
                elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_COMM
                letap[planeInc+j+edgeElems*edgeElems-edgeElems] = ghostIdx[4] +  rowInc + (j-1)
            end

            if domain.m_colLoc == 0
                elemBC[planeInc+(j-1)*edgeElems+1] |= XI_M_SYMM
            else
                elemBC[planeInc+(j-1)*edgeElems+1] |= XI_M_COMM
                lxim[planeInc+(j-1)*edgeElems+1] = ghostIdx[5] + rowInc + (j-1)
            end

            if domain.m_colLoc == domain.m_tp-1
                elemBC[planeInc+(j-1)*edgeElems+edgeElems] |= XI_P_FREE
            else
                elemBC[planeInc+(j-1)*edgeElems+edgeElems] |= XI_P_COMM
                lxip[planeInc+(j-1)*edgeElems+edgeElems] = ghostIdx[6] + rowInc + (j-1)
            end
        end
    end

    copyto!(domain.elemBC, elemBC)
    copyto!(domain.lxim, lxim)
    copyto!(domain.lxip, lxip)
    copyto!(domain.letam, letam)
    copyto!(domain.letap, letap)
    copyto!(domain.lzetam, lzetam)
    copyto!(domain.lzetap, lzetap)
end

function sortRegions(regReps::Vector{IndexT}, regSorted::Vector{IndexT}, regElemSize, numReg)
    regIndex = [v for v in 1:numReg]::Vector{IndexT}

    for i in 1:numReg-1
        for j in 1:numReg-i-1
            if regReps[j] < regReps[j+1]
                temp = regReps[j]
                regReps[j] = regReps[j+1]
                regReps[j+1] = temp

                temp = regElemSize[j]
                regElemSize[j] = regElemSize[j+1]
                regElemSize[j+1] = temp

                temp = regIndex[j]
                regIndex[j] = regIndex[j+1]
                regIndex[j+1] = temp
            end
        end
    end
    for i in 1:numReg
        regSorted[regIndex[i]] = i
    end
end

function createRegionIndexSets!(domain::Domain, nr::Int, b::Int, comm::Union{MPI.Comm, Nothing})
    domain.numReg = nr
    domain.balance = b
    @unpack_Domain domain
    myRank = getMyRank(comm)
    Random.seed!(myRank)

    regElemSize = Vector{Int}(undef, numReg)
    nextIndex::IndexT = 0

    regCSR = convert(Vector{Int}, regCSR) # records the begining and end of each region
    regReps = convert(Vector{Int}, regReps) # records the rep number per region
    regNumList = convert(Vector{IndexT}, regNumList) # Region number per domain element
    regElemlist = convert(Vector{IndexT}, regElemlist) # region indexset
    regSorted = convert(Vector{IndexT}, regSorted) # keeps index of sorted regions

    # if we only have one region just fill it
    # Fill out the regNumList with material numbers, which are always
    # the region index plus one
    if numReg == 1
        while nextIndex < numElem
            regNumList[nextIndex+1] = 1
            nextIndex+=1
        end
        regElemSize[1] = 0
    # If we have more than one region distribute the elements.
    else
        lastReg::Int = -1
        runto::IndexT = 0
        costDenominator::Int = 0
        regBinEnd = Vector{Int}(undef, numReg)
        # Determine the relative weights of all the regions.
        for i in 1:numReg
            regElemSize[i] = 0
            # INDEXING
            costDenominator += i^balance  # Total cost of all regions
            regBinEnd[i] = costDenominator  # Chance of hitting a given region is (regBinEnd[i] - regBinEdn[i-1])/costDenominator
        end
        # Until all elements are assigned
        while nextIndex < numElem
            # pick the region
            regionVar = rand(Int) % costDenominator
            i = 0
            # INDEXING
            while regionVar >= regBinEnd[i+1]
                i += 1
            end
            # rotate the regions based on MPI rank.  Rotation is Rank % NumRegions
            regionNum = ((i + myRank) % numReg) + 1
            # make sure we don't pick the same region twice in a row
            while regionNum == lastReg
                regionVar = rand(Int) % costDenominator
                i = 0
                while regionVar >= regBinEnd[i+1]
                    i += 1
                end
                regionNum = ((i + myRank) % numReg) + 1
            end
            # Pick the bin size of the region and determine the number of elements.
            binSize = rand(Int) % 1000
            if binSize < 773
                elements = rand(Int) % 15 + 1
            elseif binSize < 937
                elements = rand(Int) % 16 + 16
            elseif binSize < 970
                elements = rand(Int) % 32 + 32
            elseif binSize < 974
                elements = rand(Int) % 64 + 64
            elseif binSize < 978
                elements = rand(Int) % 128 + 128
            elseif binSize < 981
                elements = rand(Int) % 256 + 256
            else
                elements = rand(Int) % 1537 + 512
            end
            runto = elements + nextIndex
            # Store the elements.  If we hit the end before we run out of elements then just stop.
            while nextIndex < runto && nextIndex < numElem
                # INDEXING
                regNumList[nextIndex+1] = regionNum
                nextIndex += 1
            end
            lastReg = regionNum
        end
    end
    # Convert regNumList to region index sets
    # First, count size of each region
    for i in 1:numElem
        # INDEXING
        r = regNumList[i] # region index == regnum-1
        regElemSize[r]+=1
    end

    # Second, allocate each region index set
    for r in 1:numReg
        if r < div(numReg, 2)
            rep = 1
        elseif r < (numReg - div((numReg+15),20))
            rep = 1 + cost;
        else
            rep = 10 * (1 + cost)
        end
        regReps[r] = rep
    end
    sortRegions(regReps, regSorted, regElemSize, numReg);

    regCSR[1] = 0;
    # Second, allocate each region index set
    for i in 2:numReg
        regCSR[i] = regCSR[i-1] + regElemSize[i-1];
    end

    # Third, fill index sets
    for i in 1:numElem
        # INDEXING
        r = regSorted[regNumList[i]] # region index == regnum-1
        regElemlist[regCSR[r]+1] = i
        regCSR[r] += 1
    end

    # Copy to device
    copyto!(regCSR, regCSR) # records the begining and end of each region
    copyto!(regReps, regReps) # records the rep number per region
    copyto!(regNumList, regNumList) # Region number per domain element
    copyto!(regElemlist, regElemlist) # region indexset
    copyto!(regSorted, regSorted) # keeps index of sorted regions
    @pack_Domain! domain
end

function Domain(prob::LuleshProblem)
    VDF = prob.devicetype{prob.floattype}
    VDI = prob.devicetype{IndexT}
    VDInt = prob.devicetype{Int}
    colLoc = prob.col
    rowLoc = prob.row
    planeLoc = prob.plane
    nx = prob.nx
    tp = prob.side
    structured = prob.structured
    nr = prob.nr
    balance = prob.balance
    cost = prob.cost
    domain = Domain{prob.floattype}(
        prob.comm,
        VDI(), VDI(),
        VDI(), VDI(), VDI(), VDI(), VDI(), VDI(),
        VDInt(),
        VDF(), VDF(),
        VDF(),
        VDF(), VDF(), VDF(),
        VDF(),
        VDF(), VDF(), VDF(), # volo
        VDF(),
        VDF(),
        VDF(), # elemMass
        VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        VDF(), VDF(), VDF(),
        # FIXIT This is wrong
        VDF(),
        VDI(), VDI(), VDI(),
        VDInt(), VDInt(), VDI(),
        0.0, 0.0, 0.0, 0.0, 0.0, 0,
        0.0, 0.0, 0.0, 0.0, 0, 0,
        0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0,0,0,0,0,0,0,0,0,
        0,
        0,
        0,0,0,
        0,
        0,0,0, Vector{Int}(), VDInt(), VDInt(), VDI(), VDI(), VDI(),
        0,0,0,0,0,0,
        VDF(),VDF(),
        Vector{MPI.Request}(undef, 26), Vector{MPI.Request}(undef, 26)
    )

    nodelist = Vector{IndexT}()
    x = Vector{prob.floattype}()
    y = Vector{prob.floattype}()
    z = Vector{prob.floattype}()

    if structured
        domain.m_tp       = tp

        domain.m_colLoc   =   colLoc
        domain.m_rowLoc   =   rowLoc
        domain.m_planeLoc = planeLoc

        edgeElems = nx
        edgeNodes = edgeElems+1

        domain.sizeX = edgeElems
        domain.sizeY = edgeElems
        domain.sizeZ = edgeElems

        domain.numElem = domain.sizeX*domain.sizeY*domain.sizeZ ;
        # domain.padded_numElem = PAD(domain.numElem,32);

        domain.numNode = (domain.sizeX+1)*(domain.sizeY+1)*(domain.sizeZ+1)
        # domain.padded_numNode = PAD(domain.numNode,32);

        domElems = domain.numElem
        domNodes = domain.numNode
        # padded_domElems = domain.padded_numElem

        # Build domain object here. Not nice.

        allocateElemPersistent!(domain, domElems)
        allocateNodalPersistent!(domain, domNodes)

        setupCommBuffers!(domain, edgeNodes)

        initializeFields!(domain)

        buildMesh!(domain, nx, edgeNodes, edgeElems, domNodes, domElems, x, y, z, nodelist)

        domain.numSymmX = domain.numSymmY = domain.numSymmZ = 0

        if domain.m_colLoc == 0
            domain.numSymmX = (edgeElems+1)*(edgeElems+1)
        end
        if domain.m_rowLoc == 0
            domain.numSymmY = (edgeElems+1)*(edgeElems+1)
        end
        if domain.m_planeLoc == 0
            domain.numSymmZ = (edgeElems+1)*(edgeElems+1)
        end
        # Set up symmetry nodesets

        symmX = convert(Vector, domain.symmX)
        symmY = convert(Vector, domain.symmY)
        symmZ = convert(Vector, domain.symmZ)

        fill!(symmX, 0)
        fill!(symmY, 0)
        fill!(symmZ, 0)

        nidx = 1
        # INDEXING
        for i in 1:edgeNodes
            planeInc = (i-1)*edgeNodes*edgeNodes
            rowInc   = (i-1)*edgeNodes
            for j in 1:edgeNodes
                if domain.m_planeLoc == 0
                    symmZ[nidx] = rowInc   + j-1
                end
                if domain.m_rowLoc == 0
                    symmY[nidx] = planeInc + j-1
                end
                if domain.m_colLoc == 0
                    symmX[nidx] = planeInc + (j-1)*edgeNodes
                end
                nidx+=1
            end
        end
        if domain.m_planeLoc == 0
            domain.symmZ = symmZ
        end
        if domain.m_rowLoc == 0
            domain.symmY = symmY
        end
        if domain.m_colLoc == 0
            domain.symmX = symmX
        end

        setupConnectivityBC!(domain, edgeElems)
    else
        error("Reading unstructured mesh is currently missing in the Julia version of LULESH.")
    end
    # set up node-centered indexing of elements */
    nodeElemCount = zeros(IndexT, domNodes)
    # INDEXING
    for i in 1:domElems
        for j in 1:8
            nodeElemCount[nodelist[8*(i-1)+j]+1] += 1
        end
    end

    nodeElemStart = zeros(IndexT, domNodes)
    nodeElemStart[1] = 0
    for i in 2:domNodes
        nodeElemStart[i] = nodeElemStart[i-1] + nodeElemCount[i-1]
    end
    nodeElemCornerList = Vector{IndexT}(undef, nodeElemStart[domNodes] + nodeElemCount[domNodes] )

    nodeElemCount .= 0

    for i in 1:domElems
        for j in 1:8
            # @show i,j
            m = nodelist[8*(i-1)+j] + 1
            k = 8*(i-1) + j
            # INDEXING
            offset = nodeElemStart[m] + nodeElemCount[m]
            nodeElemCornerList[offset+1] = k
            nodeElemCount[m] += 1
        end
    end

    clSize = nodeElemStart[domNodes] + nodeElemCount[domNodes]
    for i in 1:clSize
        clv = nodeElemCornerList[i] ;
        if (clv < 0) || (clv > domElems*8)
            error("AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!")
        end
    end

    domain.nodeElemStart = convert(VDI, nodeElemStart)
    domain.nodeElemCount = convert(VDI, nodeElemCount)
    domain.nodeElemCornerList = convert(VDI, nodeElemCornerList)

    # Create a material IndexSet (entire domain same material for now)
    matElemlist = Vector{IndexT}(undef, domElems)
    for i in 1:domElems
        matElemlist[i] = i
    end
    copyto!(domain.matElemlist, matElemlist)

    domain.bad_vol = -1
    domain.bad_q = -1
    domain.dthydro = 1e20
    domain.dtcourant = 1e20

    # initialize material parameters
    domain.time      = 0.
    domain.dtfixed = -1.0e-6
    domain.deltatimemultlb = 1.1
    domain.deltatimemultub = 1.2
    domain.stoptime  = 1.0e-2
    domain.dtmax     = 1.0e-2
    domain.cycle   = 0

    domain.e_cut = 1.0e-7
    domain.p_cut = 1.0e-7
    domain.q_cut = 1.0e-7
    domain.u_cut = 1.0e-7
    domain.v_cut = 1.0e-10

    domain.hgcoef      = 3.0
    domain.ss4o3       = 4.0/3.0

    domain.qstop              =  1.0e+12
    domain.monoq_max_slope    =  1.0
    domain.monoq_limiter_mult =  2.0
    domain.qlc_monoq          = 0.5
    domain.qqc_monoq          = 2.0/3.0
    domain.qqc                = 2.0

    domain.pmin =  0.
    domain.emin = -1.0e+15

    domain.dvovmax =  0.1

    domain.eosvmax =  1.0e+9
    domain.eosvmin =  1.0e-9

    domain.refdens =  1.0

    # initialize field data
    nodalMass = Vector{prob.floattype}(undef, domNodes)
    volo = Vector{prob.floattype}(undef, domElems)
    elemMass = Vector{prob.floattype}(undef, domElems)
    fill!(nodalMass, 0)

    for i in 1:domElems
        x_local = Vector{prob.floattype}(undef, 8)
        y_local = Vector{prob.floattype}(undef, 8)
        z_local = Vector{prob.floattype}(undef, 8)
        for lnode in 1:8
            gnode = nodelist[(i-1)*8+lnode]+1
            x_local[lnode] = x[gnode]
            y_local[lnode] = y[gnode]
            z_local[lnode] = z[gnode]
        end
        # volume calculations
        volume = calcElemVolume(x_local, y_local, z_local )
        volo[i] = volume
        elemMass[i] = volume
        for j in 1:8
            gnode = nodelist[(i-1)*8+j]+1
            nodalMass[gnode] += volume / 8.0
        end
    end

    copyto!(domain.nodalMass, nodalMass)
    copyto!(domain.volo, volo)
    copyto!(domain.elemMass, elemMass)

    # deposit energy
    domain.octantCorner = 0;
    # deposit initial energy
    # An energy of 3.948746e+7 is correct for a problem with
    # 45 zones along a side - we need to scale it
    ebase = 3.948746e+7
    scale = (nx*domain.m_tp)/45.0;
    einit = ebase*scale*scale*scale;
    if domain.m_rowLoc + domain.m_colLoc + domain.m_planeLoc == 0
        # Dump into the first zone (which we know is in the corner)
        # of the domain that sits at the origin
        # TODO This only works for CUDA
        domain.e[1] = einit;
    end

    # set initial deltatime base on analytic CFL calculation
    domain.deltatime = (.5*cbrt(domain.volo[1]))/sqrt(2.0*einit)

    domain.cost = cost
    resize!(domain.regNumList, domain.numElem)  # material indexset
    resize!(domain.regElemlist, domain.numElem)  # material indexset
    resize!(domain.regCSR, nr)
    resize!(domain.regReps, nr)
    resize!(domain.regSorted, nr)

    # Setup region index sets. For now, these are constant sized
    # throughout the run, but could be changed every cycle to
    # simulate effects of ALE on the lagrange solver
    createRegionIndexSets!(domain, nr, balance, prob.comm)
    # Setup symmetry nodesets
    setupSymmetryPlanes(domain, edgeNodes)

    # Setup element connectivities
    setupElementConnectivities(domain, edgeElems)

    # Setup symmetry planes and free surface boundary arrays
    setupBoundaryConditions(domain, edgeElems)
    return domain
end

function  setupSymmetryPlanes(domain, edgeNodes)
    nidx = 1
    for i in 0:(edgeNodes-1)
        planeInc = i*edgeNodes*edgeNodes
        rowInc   = i*edgeNodes
        for j in 0:(edgeNodes-1)
            if domain.m_planeLoc == 0
                domain.symmZ[nidx] = rowInc   + j
            end
            if domain.m_rowLoc == 0
                domain.symmY[nidx] = planeInc + j
            end
            if domain.m_colLoc == 0
                domain.symmX[nidx] = planeInc + j*edgeNodes
            end
            nidx += 1
        end
    end
end

function setupElementConnectivities(domain::Domain, edgeElems)
    domain.lxim[1] = 1
    for i in 1:(domain.numElem - 1)
        domain.lxim[i+1] = i-1
        domain.lxip[i] = i
    end
    domain.lxip[domain.numElem] = domain.numElem - 1

    for i in 0:(edgeElems - 1)
        domain.letam[i+1] = i
        domain.letap[domain.numElem - edgeElems + i + 1] = i
    end

    for i in edgeElems:(domain.numElem - 1)
        domain.letam[i+1] = i - edgeElems
        domain.letap[i-edgeElems + 1] = i
    end

    for i in 0:(edgeElems*edgeElems - 1)
        domain.lzetam[i+1] = i
        domain.lzetap[domain.numElem - edgeElems*edgeElems + i + 1] = domain.numElem - edgeElems*edgeElems+i
    end

    for i in edgeElems*edgeElems:(domain.numElem - 1)
        domain.lzetam[i+1] = i - edgeElems * edgeElems
        domain.lzetap[i - edgeElems*edgeElems+1] = i
    end
end

function setupBoundaryConditions(domain::Domain, edgeElems)
  ghostIdx = Vector{IndexT}(undef, 6) # offsets to ghost locations

  # set up boundary condition information
  for i in 0:domain.numElem - 1
    domain.elemBC[i+1] = 0
  end

  for i in 1:6
    ghostIdx[i] = typemin(IndexT)
  end

  pidx = domain.numElem

  if domain.m_planeMin != 0
    ghostIdx[1] = pidx
    pidx += domain.sizeX*domain.sizeY
  end

  if m_planeMax != 0
    ghostIdx[2] = pidx
    pidx += domain.sizeX*domain.sizeY
  end

  if m_rowMin != 0
    ghostIdx[3] = pidx
    pidx += domain.sizeX*domain.sizeZ
  end

  if m_rowMax != 0
    ghostIdx[4] = pidx
    pidx += domain.sizeX*domain.sizeZ
  end

  if m_colMin != 0
    ghostIdx[5] = pidx
    pidx += domain.sizeY*domain.sizeZ
  end

  if m_colMax != 0
    ghostIdx[6] = pidx
  end


  # symmetry plane or free surface BCs

    for i in 1:edgeElems
        planeInc = (i-1)*edgeElems*edgeElems
        rowInc   = (i-1)*edgeElems
        for j in 1:edgeElems
            if domain.m_planeLoc == 0
                domain.elemBC[rowInc+j] |= ZETA_M_SYMM
            else
                domain.elemBC[rowInc+j] |= ZETA_M_COMM
                domain.lzetam[rowInc+j] = ghostIdx[1] + rowInc + (j-1)
            end

            if domain.m_planeLoc == domain.m_tp-1
                domain.elemBC[rowInc+j+domain.numElem-edgeElems*edgeElems] |= ZETA_P_FREE
            else
                domain.elemBC[rowInc+j+domain.numElem-edgeElems*edgeElems] |= ZETA_P_COMM
                domain.lzetap[rowInc+j+domain.numElem-edgeElems*edgeElems] = ghostIdx[2] + rowInc + (j-1)
            end

            if domain.m_rowLoc == 0
                domain.elemBC[planeInc+j] |= ETA_M_SYMM
            else
                domain.elemBC[planeInc+j] |= ETA_M_COMM
                domain.letam[planeInc+j] = ghostIdx[3] + rowInc + (j-1)
            end


            if domain.m_rowLoc == domain.m_tp-1
                domain.elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_FREE
            else
                domain.elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_COMM
                domain.letap[planeInc+j+edgeElems*edgeElems-edgeElems] = ghostIdx[4] +  rowInc + (j-1)
            end

            if domain.m_colLoc == 0
                domain.elemBC[planeInc+(j-1)*edgeElems+1] |= XI_M_SYMM
            else
                domain.elemBC[planeInc+(j-1)*edgeElems+1] |= XI_M_COMM
                domain.lxim[planeInc+(j-1)*edgeElems+1] = ghostIdx[5] + rowInc + (j-1)
            end

            if domain.m_colLoc == domain.m_tp-1
                domain.elemBC[planeInc+(j-1)*edgeElems+edgeElems] |= XI_P_FREE
            else
                domain.elemBC[planeInc+(j-1)*edgeElems+edgeElems] |= XI_P_COMM
                domain.lxip[planeInc+(j-1)*edgeElems+edgeElems] = ghostIdx[6] + rowInc + (j-1)
            end
        end
    end
end



function timeIncrement!(domain::Domain)
    targetdt = domain.stoptime - domain.time
    if domain.dtfixed <= 0.0 && domain.cycle != 0
        olddt = domain.deltatime

        gnewdt = typemax(Float64)
        newdt = 0.0

        if domain.dtcourant < gnewdt
            gnewdt = domain.dtcourant / 2.0
        end

        if domain.dthydro < gnewdt
            gnewdt = domain.dthydro * 2.0 / 3.0
        end
        newdt = comm_min(gnewdt, domain.comm)

        ratio = newdt / olddt
        if ratio >= 1.0
            if ratio < domain.deltatimemultlb
                newdt = olddt
            elseif ratio > domain.deltatimemultub
                newdt = olddt * domain.deltatimemultub
            end
        end

        newdt = min(newdt, domain.dtmax)

        domain.deltatime = newdt
    end

    # try to prevent very small scaling on the next cycle
    if domain.deltatime < targetdt < 4.0 * domain.deltatime / 3.0
        targetdt = 2.0 * domain.deltatime / 3.0
    end

    if targetdt < domain.deltatime
        domain.deltatime = targetdt
    end
    domain.time += domain.deltatime
    domain.cycle += 1
end

function initStressTermsForElems(domain::Domain, sigxx, sigyy, sigzz)
    # Based on FORTRAN implementation down from here
    @assert axes(sigxx) == axes(sigyy) == axes(sigzz) == axes(domain.p) == axes(domain.q)
    for i in 1:domain.numElem
        sigxx[i] = sigyy[i] = sigzz[i] = - domain.p[i] - domain.q[i]
    end
end

function sumElemFaceNormal(x0,  y0,  z0,
                           x1,  y1,  z1,
                           x2,  y2,  z2,
                           x3,  y3,  z3)
  RHALF = 0.5
  RQTR = 0.25

  bisectX0 = RHALF * (x3 + x2 - x1 - x0)
  bisectY0 = RHALF * (y3 + y2 - y1 - y0)
  bisectZ0 = RHALF * (z3 + z2 - z1 - z0)
  bisectX1 = RHALF * (x2 + x1 - x3 - x0)
  bisectY1 = RHALF * (y2 + y1 - y3 - y0)
  bisectZ1 = RHALF * (z2 + z1 - z3 - z0)
  areaX = RQTR * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1)
  areaY = RQTR * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1)
  areaZ = RQTR * (bisectX0 * bisectY1 - bisectY0 * bisectX1)

  return areaX, areaY, areaZ
end

@inline function calcElemNodeNormals(x, y, z)
    @inbounds @views begin
    pf = zeros(MMatrix{8, 3, Float64})
    pfx = view(pf, :, 1)
    pfy = view(pf, :, 2)
    pfz = view(pf, :, 3)

    # evaluate face one: nodes 1, 2, 3, 4
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[1], y[1], z[1], x[2], y[2], z[2],
                            x[3], y[3], z[3], x[4], y[4], z[4])
    pfx[1:4] .+= areaX
    pfy[1:4] .+= areaY
    pfz[1:4] .+= areaZ

    # evaluate face two: nodes 1, 5, 6, 2
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[1], y[1], z[1], x[5], y[5], z[5],
                            x[6], y[6], z[6], x[2], y[2], z[2])

    pfx[1:2] .+= areaX
    pfx[5:6] .+= areaX
    pfy[1:2] .+= areaY
    pfy[5:6] .+= areaY
    pfz[1:2] .+= areaZ
    pfz[5:6] .+= areaZ

    #evaluate face three: nodes 2, 6, 7, 3
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[2], y[2], z[2], x[6], y[6], z[6],
                            x[7], y[7], z[7], x[3], y[3], z[3])

    pfx[2:3] .+= areaX
    pfx[6:7] .+= areaX
    pfy[2:3] .+= areaY
    pfy[6:7] .+= areaY
    pfz[2:3] .+= areaZ
    pfz[6:7] .+= areaZ

    #evaluate face four: nodes 3, 7, 8, 4
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[3], y[3], z[3], x[7], y[7], z[7],
                            x[8], y[8], z[8], x[4], y[4], z[4])

    pfx[3:4] .+= areaX
    pfx[7:8] .+= areaX
    pfy[3:4] .+= areaY
    pfy[7:8] .+= areaY
    pfz[3:4] .+= areaZ
    pfz[7:8] .+= areaZ

    # evaluate face five: nodes 4, 8, 5, 1
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[4], y[4], z[4], x[8], y[8], z[8],
                            x[5], y[5], z[5], x[1], y[1], z[1])

    pfx[1]    += areaX
    pfx[4:5] .+= areaX
    pfx[8]    += areaX
    pfy[1]    += areaY
    pfy[4:5] .+= areaY
    pfy[8]    += areaY
    pfz[1]    += areaZ
    pfz[4:5] .+= areaZ
    pfz[8]    += areaZ

    # evaluate face six: nodes 5, 8, 7, 6
    areaX, areaY, areaZ = sumElemFaceNormal(
                            x[5], y[5], z[5], x[8], y[8], z[8],
                            x[7], y[7], z[7], x[6], y[6], z[6])
    pfx[5:8] .+= areaX
    pfy[5:8] .+= areaY
    pfz[5:8] .+= areaZ

    return SMatrix(pf)
    end #@inbounds
end

@inline function sumElemStressesToNodeForces(B, sig_xx, sig_yy, sig_zz,  fx_elem,  fy_elem,  fz_elem, k)

  @inbounds begin
    stress_xx = sig_xx[k]
    stress_yy = sig_yy[k]
    stress_zz = sig_zz[k]

    fx = -stress_xx .* B[:, 1]
    fy = -stress_yy .* B[:, 2]
    fz = -stress_zz .* B[:, 3]

    fx_elem[(k-1)*8+1:k*8] = fx
    fy_elem[(k-1)*8+1:k*8] = fy
    fz_elem[(k-1)*8+1:k*8] = fz
  end
end

function JcollectNodal(i, nodelist, sources, src_x)

    src_x_shared = JACC.shared(src_x)
    # myARRAY = zeros(Float64, 8)
    @inbounds begin

        i8 = (i - 1) * 8
        sources[i8 + 1] = src_x[nodelist[i8 + 1] + 1]
        sources[i8 + 2] = src_x[nodelist[i8 + 2] + 1]
        sources[i8 + 3] = src_x[nodelist[i8 + 3] + 1]
        sources[i8 + 4] = src_x[nodelist[i8 + 4] + 1]
        sources[i8 + 5] = src_x[nodelist[i8 + 5] + 1]
        sources[i8 + 6] = src_x[nodelist[i8 + 6] + 1]
        sources[i8 + 7] = src_x[nodelist[i8 + 7] + 1]
        sources[i8 + 8] = src_x[nodelist[i8 + 8] + 1]

    end

end

function JcalcElemShapeFunctionDerivatives(i, x, y, z, b, determ)
    @inbounds begin

        i8 = (i - 1) * 8
        i3 = (i - 1) * 24


        x0 = x[i8 + 1]
        x1 = x[i8 + 2]
        x2 = x[i8 + 3]
        x3 = x[i8 + 4]
        x4 = x[i8 + 5]
        x5 = x[i8 + 6]
        x6 = x[i8 + 7]
        x7 = x[i8 + 8]

        y0 = y[i8 + 1]
        y1 = y[i8 + 2]
        y2 = y[i8 + 3]
        y3 = y[i8 + 4]
        y4 = y[i8 + 5]
        y5 = y[i8 + 6]
        y6 = y[i8 + 7]
        y7 = y[i8 + 8]

        z0 = z[i8 + 1]
        z1 = z[i8 + 2]
        z2 = z[i8 + 3]
        z3 = z[i8 + 4]
        z4 = z[i8 + 5]
        z5 = z[i8 + 6]
        z6 = z[i8 + 7]
        z7 = z[i8 + 8]



        fjxxi = .125 * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) )
        fjxet = .125 * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) )
        fjxze = .125 * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) )

        fjyxi = .125 * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) )
        fjyet = .125 * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) )
        fjyze = .125 * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) )

        fjzxi = .125 * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) )
        fjzet = .125 * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) )
        fjzze = .125 * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) )

        # compute cofactors
        cjxxi =    (fjyet * fjzze) - (fjzet * fjyze)
        cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze)
        cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet)

        cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze)
        cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze)
        cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet)

        cjzxi =    (fjxet * fjyze) - (fjyet * fjxze)
        cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze)
        cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet)

        # calculate partials :
        #     this need only be done for l = 0,1,2,3   since , by symmetry ,
        #     (6,7,4,5) = - (0,1,2,3) .
        b[1 + i3] =   -  cjxxi  -  cjxet  -  cjxze
        b[2 + i3] =      cjxxi  -  cjxet  -  cjxze
        b[3 + i3] =      cjxxi  +  cjxet  -  cjxze
        b[4 + i3] =   -  cjxxi  +  cjxet  -  cjxze
        b[5 + i3] = -b[3 + i3]
        b[6 + i3] = -b[4 + i3]
        b[7 + i3] = -b[1 + i3]
        b[8 + i3] = -b[2 + i3]

        b[9 + i3] =   -  cjyxi  -  cjyet  -  cjyze
        b[10 + i3] =      cjyxi  -  cjyet  -  cjyze
        b[11 + i3] =      cjyxi  +  cjyet  -  cjyze
        b[12 + i3] =   -  cjyxi  +  cjyet  -  cjyze
        b[13 + i3] = -b[11 + i3]
        b[14 + i3] = -b[12 + i3]
        b[15 + i3] = -b[9 + i3]
        b[16 + i3] = -b[10 + i3]

        b[17 + i3] =   -  cjzxi  -  cjzet  -  cjzze
        b[18 + i3] =      cjzxi  -  cjzet  -  cjzze
        b[19 + i3] =      cjzxi  +  cjzet  -  cjzze
        b[20 + i3] =   -  cjzxi  +  cjzet  -  cjzze
        b[21 + i3] = -b[19 + i3]
        b[22 + i3] = -b[20 + i3]
        b[23 + i3] = -b[17 + i3]
        b[24 + i3] = -b[18 + i3]


    end #inbounds

    # calculate jacobian determinant (volume)
    determ[i] = 8.0 * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet)
end



function JupdatedcalcElemNodeNormals(i, x, y, z, pf)
    @inbounds begin
        # if i > 27000
        #     return
        # end
        i3 = (i - 1) * 8
        i8 = (i - 1) * 24


        # Helper function to update pf values
        function update_pf!(offset, areaX, areaY, areaZ)
            pf[offset + i8] += areaX
            pf[offset + i8 + 8] += areaY
            pf[offset + i8 + 16] += areaZ
        end
        # Evaluate face one: nodes 1, 2, 3, 4
        areaX, areaY, areaZ = sumElemFaceNormal(
            x[1+i3], y[1+i3], z[1+i3], x[2+i3], y[2+i3], z[2+i3],
            x[3+i3], y[3+i3], z[3+i3], x[4+i3], y[4+i3], z[4+i3])
        
        for j in 1:4
            update_pf!(j, areaX, areaY, areaZ)
        end

        # Evaluate face two: nodes 1, 5, 6, 2
        areaX, areaY, areaZ = sumElemFaceNormal(
            x[1+i3], y[1+i3], z[1+i3], x[5+i3], y[5+i3], z[5+i3],
            x[6+i3], y[6+i3], z[6+i3], x[2+i3], y[2+i3], z[2+i3])
        
        update_pf!(1, areaX, areaY, areaZ)
        update_pf!(2, areaX, areaY, areaZ)
        update_pf!(5, areaX, areaY, areaZ)
        update_pf!(6, areaX, areaY, areaZ)

        # Evaluate face three: nodes 2, 6, 7, 3
        areaX, areaY, areaZ = sumElemFaceNormal(
            x[2+i3], y[2+i3], z[2+i3], x[6+i3], y[6+i3], z[6+i3],
            x[7+i3], y[7+i3], z[7+i3], x[3+i3], y[3+i3], z[3+i3])
        
        update_pf!(2, areaX, areaY, areaZ)
        update_pf!(3, areaX, areaY, areaZ)
        update_pf!(6, areaX, areaY, areaZ)
        update_pf!(7, areaX, areaY, areaZ)

        # Evaluate face four: nodes 3, 7, 8, 4
        areaX, areaY, areaZ = sumElemFaceNormal(
            x[3+i3], y[3+i3], z[3+i3], x[7+i3], y[7+i3], z[7+i3],
            x[8+i3], y[8+i3], z[8+i3], x[4+i3], y[4+i3], z[4+i3])
        
        update_pf!(3, areaX, areaY, areaZ)
        update_pf!(4, areaX, areaY, areaZ)
        update_pf!(7, areaX, areaY, areaZ)
        update_pf!(8, areaX, areaY, areaZ)

        # Evaluate face five: nodes 4, 8, 5, 1
        areaX, areaY, areaZ = sumElemFaceNormal(
            x[4+i3], y[4+i3], z[4+i3], x[8+i3], y[8+i3], z[8+i3],
            x[5+i3], y[5+i3], z[5+i3], x[1+i3], y[1+i3], z[1+i3])
        
        update_pf!(1, areaX, areaY, areaZ)
        update_pf!(4, areaX, areaY, areaZ)
        update_pf!(5, areaX, areaY, areaZ)
        update_pf!(8, areaX, areaY, areaZ)

        # Evaluate face six: nodes 5, 8, 7, 6
        areaX, areaY, areaZ = sumElemFaceNormal(
            x[5+i3], y[5+i3], z[5+i3], x[8+i3], y[8+i3], z[8+i3],
            x[7+i3], y[7+i3], z[7+i3], x[6+i3], y[6+i3], z[6+i3])
        
        for j in 5:8
            update_pf!(j, areaX, areaY, areaZ)
        end
    end
end

function JsumElemStressesToNodeForces(i, B, sig_xx, sig_yy, sig_zz,  fx_elem,  fy_elem,  fz_elem, fx, fy, fz)

    @inbounds begin

        i3 = (i - 1) * 8
        i8 = (i - 1) * 24
        stress_xx = sig_xx[i]
        stress_yy = sig_yy[i]
        stress_zz = sig_zz[i]

        for j in 1:8
            fx[j + i3] = -stress_xx * B[j + (i-1)*24]  # Assuming B was originally 8x3xnumElem
            fy[j + i3] = -stress_yy * B[j + 8 + (i-1)*24]
            fz[j + i3] = -stress_zz * B[j + 16 + (i-1)*24]
        end

        for j in 1:8
            fx_elem[i3 + j] = fx[j + i3]
            fy_elem[i3 + j] = fy[j + i3]
            fz_elem[i3 + j] = fz[j + i3]
        end
    end
end


function Jupdate_domain(i, fx_arr, fy_arr, fz_arr, nodeElemCount, nodeElemStart, nodeElemCornerList, fx_elem, fy_elem, fz_elem)

    @inbounds begin
        # if i > 29791
        #     return
        # end
        count = nodeElemCount[i]
        start = nodeElemStart[i]
        fx = 0.0
        fy = 0.0
        fz = 0.0
        for j in 1:count
            elem = nodeElemCornerList[start + j]
            fx += fx_elem[elem]
            fy += fy_elem[elem]
            fz += fz_elem[elem]
        end

        fx_arr[i] += fx
        fy_arr[i] += fy
        fz_arr[i] += fz

    end

end

function MultiDintegrateStressForElems(domain::Domain, sigxx, sigyy, sigzz, determ)
    # Based on FORTRAN implementation down from here
    # loop over all elements
    numElem8 = domain.numElem*8
    numElem3 = domain.numElem*3
    T = typeof(domain.x)
    
    fx_elem = T(undef, numElem8)
    fy_elem = T(undef, numElem8)
    fz_elem = T(undef, numElem8)
    # FIXIT. This has to be device type
    nodelist = domain.nodelist

    # One Dimensional offsets
    x_single = Vector{Float64}(undef, numElem8)
    y_single = Vector{Float64}(undef, numElem8)
    z_single = Vector{Float64}(undef, numElem8)
    b_single = Vector{Float64}(undef, numElem3 * 8)
    pf_single = Vector{Float64}(undef, numElem8 * 3)
    pfx_single = Vector{Float64}(undef, numElem8)
    pfy_single = Vector{Float64}(undef, numElem8)
    pfz_single = Vector{Float64}(undef, numElem8)
    fx_single = Vector{Float64}(undef, numElem8)
    fy_single = Vector{Float64}(undef, numElem8)
    fz_single = Vector{Float64}(undef, numElem8)

    x_domain_d = JACC.Array(domain.x)
    y_domain_d = JACC.Array(domain.y)
    z_domain_d = JACC.Array(domain.z)
    nodelist_d = JACC.Array(nodelist)
    determ_d = JACC.Array(determ)
    sig_xx_d = JACC.Array(sigxx)
    sig_yy_d = JACC.Array(sigyy)
    sig_zz_d = JACC.Array(sigzz)
    fx_elem_d = JACC.Array(fx_elem)
    fy_elem_d = JACC.Array(fy_elem)
    fz_elem_d = JACC.Array(fz_elem)
    fx_d = JACC.Array(domain.fx)
    fy_d = JACC.Array(domain.fy)
    fz_d = JACC.Array(domain.fz)
    nodeElemCount_d = JACC.Array(domain.nodeElemCount)
    nodeElemStart_d = JACC.Array(domain.nodeElemStart)
    nodeElemCornerList_d = JACC.Array(domain.nodeElemCornerList)

    x_single_d = JACC.Array(x_single)
    y_single_d = JACC.Array(y_single)
    z_single_d = JACC.Array(z_single)
    fx_single_d = JACC.Array(fx_single)
    fy_single_d = JACC.Array(fy_single)
    fz_single_d = JACC.Array(fz_single)
    b_single_d = JACC.Array(b_single)
    pf_single_d = JACC.Array(pf_single)

    numNode = domain.numNode

    time1 = time()
    
    JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, x_single_d, x_domain_d)
    
    JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, y_single_d, y_domain_d)
    
    JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, z_single_d, z_domain_d)
    
    JACC.parallel_for(domain.numElem, JcalcElemShapeFunctionDerivatives, x_single_d, y_single_d, z_single_d, b_single_d, determ_d)

    
    
    JACC.parallel_for(domain.numElem, JupdatedcalcElemNodeNormals, x_single_d, y_single_d, z_single_d, pf_single_d)
    
    JACC.parallel_for(domain.numElem, JsumElemStressesToNodeForces, pf_single_d, sig_xx_d, sig_yy_d, sig_zz_d, fx_elem_d, fy_elem_d, fz_elem_d, fx_single_d, fy_single_d, fz_single_d)
    
    JACC.parallel_for(numNode, Jupdate_domain, fx_d, fy_d, fz_d, nodeElemCount_d, nodeElemStart_d, nodeElemCornerList_d, fx_elem_d, fy_elem_d, fz_elem_d)
    
    time2 = time()

    println("Time taken for warmup run integrateStressForElems: ", time2 - time1)

    total_time_collectNodal = 0.0
    total_time_calcElemShapeFunctionDerivatives = 0.0
    total_time_updatedcalcElemNodeNormals = 0.0
    total_time_sumElemStressesToNodeForces = 0.0
    total_time_update_domain = 0.0
    num_runs = 10

    for i in 1:num_runs
        time1 = time()
        JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, x_single_d, x_domain_d)
        JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, y_single_d, y_domain_d)
        JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, z_single_d, z_domain_d)
        time2 = time()
        total_time_collectNodal += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, JcalcElemShapeFunctionDerivatives, x_single_d, y_single_d, z_single_d, b_single_d, determ_d)
        time2 = time()
        total_time_calcElemShapeFunctionDerivatives += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, JupdatedcalcElemNodeNormals, x_single_d, y_single_d, z_single_d, pf_single_d)
        time2 = time()
        total_time_updatedcalcElemNodeNormals += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, JsumElemStressesToNodeForces, pf_single_d, sig_xx_d, sig_yy_d, sig_zz_d, fx_elem_d, fy_elem_d, fz_elem_d, fx_single_d, fy_single_d, fz_single_d)
        time2 = time()
        total_time_sumElemStressesToNodeForces += (time2 - time1)

        time1 = time()
        JACC.parallel_for(numNode, Jupdate_domain, fx_d, fy_d, fz_d, nodeElemCount_d, nodeElemStart_d, nodeElemCornerList_d, fx_elem_d, fy_elem_d, fz_elem_d)
        time2 = time()
        total_time_update_domain += (time2 - time1)
    end

    average_time_collectNodal = total_time_collectNodal / num_runs
    average_time_calcElemShapeFunctionDerivatives = total_time_calcElemShapeFunctionDerivatives / num_runs
    average_time_updatedcalcElemNodeNormals = total_time_updatedcalcElemNodeNormals / num_runs
    average_time_sumElemStressesToNodeForces = total_time_sumElemStressesToNodeForces / num_runs
    average_time_update_domain = total_time_update_domain / num_runs

    total_average_time = average_time_collectNodal + average_time_calcElemShapeFunctionDerivatives + average_time_updatedcalcElemNodeNormals + average_time_sumElemStressesToNodeForces + average_time_update_domain

    println("Average time taken for JcollectNodal over $num_runs runs: ", average_time_collectNodal)
    println("Average time taken for JcalcElemShapeFunctionDerivatives over $num_runs runs: ", average_time_calcElemShapeFunctionDerivatives)
    println("Average time taken for JupdatedcalcElemNodeNormals over $num_runs runs: ", average_time_updatedcalcElemNodeNormals)
    println("Average time taken for JsumElemStressesToNodeForces over $num_runs runs: ", average_time_sumElemStressesToNodeForces)
    println("Average time taken for Jupdate_domain over $num_runs runs: ", average_time_update_domain)
    println("Total average time taken for all functions over $num_runs runs: ", total_average_time)

    copyto!(determ, determ_d)
    for i in 1:domain.numElem
        if determ[i] <= 0.0
            println("Value and index of error: ", determ[i], ": ", i)
            error("Early Volume Error")
        end
    end

    # # Use these with JACC
    copyto!(domain.fx, fx_d)
    copyto!(domain.fy, fy_d)
    copyto!(domain.fz, fz_d)



end


function integrateStressForElems(domain::Domain, sigxx, sigyy, sigzz, determ)
    # Based on FORTRAN implementation down from here
    # loop over all elements
    numElem8 = domain.numElem*8
    T = typeof(domain.x)
    fx_elem = T(undef, numElem8)
    fy_elem = T(undef, numElem8)
    fz_elem = T(undef, numElem8)
    # FIXIT. This has to be device type
    nodelist = domain.nodelist
    @inbounds for k in 1:domain.numElem
        x_local = collectNodal(nodelist, domain.x, (k-1)*8)
        y_local = collectNodal(nodelist, domain.y, (k-1)*8)
        z_local = collectNodal(nodelist, domain.z, (k-1)*8)
        _, detJ = calcElemShapeFunctionDerivatives(x_local, y_local, z_local)
        determ[k] = detJ
        if determ[k] <= 0.0
            error("Early Volume Error")
        end
        B = calcElemNodeNormals(x_local, y_local, z_local)
        sumElemStressesToNodeForces(B, sigxx, sigyy, sigzz, fx_elem, fy_elem, fz_elem, k)
    end

    numNode = domain.numNode

    @inbounds for gnode in 1:numNode
        count = domain.nodeElemCount[gnode]
        start = domain.nodeElemStart[gnode]
        fx = zero(eltype(fx_elem))
        fy = zero(eltype(fy_elem))
        fz = zero(eltype(fz_elem))
        @simd for i in 1:count
            elem = domain.nodeElemCornerList[start+i]
            fx = fx + fx_elem[elem]
            fy = fy + fy_elem[elem]
            fz = fz + fz_elem[elem]
        end
        domain.fx[gnode] = fx
        domain.fy[gnode] = fy
        domain.fz[gnode] = fz
    end
end

@inline function collectDomainNodesToElemNodes(domain::Domain, i)
    @inbounds begin
    nd0i = domain.nodelist[i]
    nd1i = domain.nodelist[i+1]
    nd2i = domain.nodelist[i+2]
    nd3i = domain.nodelist[i+3]
    nd4i = domain.nodelist[i+4]
    nd5i = domain.nodelist[i+5]
    nd6i = domain.nodelist[i+6]
    nd7i = domain.nodelist[i+7]

    elemX = SVector(
        domain.x[nd0i+1],
        domain.x[nd1i+1],
        domain.x[nd2i+1],
        domain.x[nd3i+1],
        domain.x[nd4i+1],
        domain.x[nd5i+1],
        domain.x[nd6i+1],
        domain.x[nd7i+1],
    )

    elemY = SVector(
        domain.y[nd0i+1],
        domain.y[nd1i+1],
        domain.y[nd2i+1],
        domain.y[nd3i+1],
        domain.y[nd4i+1],
        domain.y[nd5i+1],
        domain.y[nd6i+1],
        domain.y[nd7i+1],
    )

    elemZ = SVector(
       domain.z[nd0i+1],
       domain.z[nd1i+1],
       domain.z[nd2i+1],
       domain.z[nd3i+1],
       domain.z[nd4i+1],
       domain.z[nd5i+1],
       domain.z[nd6i+1],
       domain.z[nd7i+1],
    )

    return elemX, elemY, elemZ
    end
end

@inline function voluDer(x0, x1, x2,
                   x3, x4, x5,
                   y0, y1, y2,
                   y3, y4, y5,
                   z0, z1, z2,
                   z3, z4, z5)

    dvdx =
        (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +
        (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -
        (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5)

    dvdy =
        - (x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) -
        (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +
        (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5)

    dvdz =
        - (y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) -
        (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +
        (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5)

    twelfth = 1.0 / 12.0

    dvdx = dvdx * twelfth
    dvdy = dvdy * twelfth
    dvdz = dvdz * twelfth

    return dvdx, dvdy, dvdz
end

@inline function calcElemVolumeDerivative(x, y, z)

    dvdx1, dvdy1, dvdz1 = voluDer(
        x[2], x[3], x[4], x[5], x[6], x[8],
        y[2], y[3], y[4], y[5], y[6], y[8],
        z[2], z[3], z[4], z[5], z[6], z[8])

    dvdx2, dvdy2, dvdz2 = voluDer(
        x[1], x[2], x[3], x[8], x[5], x[7],
        y[1], y[2], y[3], y[8], y[5], y[7],
        z[1], z[2], z[3], z[8], z[5], z[7])

    dvdx3, dvdy3, dvdz3 = voluDer(
        x[4], x[1], x[2], x[7], x[8], x[6],
        y[4], y[1], y[2], y[7], y[8], y[6],
        z[4], z[1], z[2], z[7], z[8], z[6])

    dvdx4, dvdy4, dvdz4 = voluDer(
        x[3], x[4], x[1], x[6], x[7], x[5],
        y[3], y[4], y[1], y[6], y[7], y[5],
        z[3], z[4], z[1], z[6], z[7], z[5])

    dvdx5, dvdy5, dvdz5 = voluDer(
        x[8], x[7], x[6], x[1], x[4], x[2],
        y[8], y[7], y[6], y[1], y[4], y[2],
        z[8], z[7], z[6], z[1], z[4], z[2])

    dvdx6, dvdy6, dvdz6 = voluDer(
        x[5], x[8], x[7], x[2], x[1], x[3],
        y[5], y[8], y[7], y[2], y[1], y[3],
        z[5], z[8], z[7], z[2], z[1], z[3])

    dvdx7, dvdy7, dvdz7 = voluDer(
        x[6], x[5], x[8], x[3], x[2], x[4],
        y[6], y[5], y[8], y[3], y[2], y[4],
        z[6], z[5], z[8], z[3], z[2], z[4])

    dvdx8, dvdy8, dvdz8 = voluDer(
        x[7], x[6], x[5], x[4], x[3], x[1],
        y[7], y[6], y[5], y[4], y[3], y[1],
        z[7], z[6], z[5], z[4], z[3], z[1])

    dvdx = SVector(dvdx1, dvdx2, dvdx3, dvdx4, dvdx5, dvdx6, dvdx7, dvdx8)
    dvdy = SVector(dvdy1, dvdy2, dvdy3, dvdy4, dvdy5, dvdy6, dvdy7, dvdy8)
    dvdz = SVector(dvdz1, dvdz2, dvdz3, dvdz4, dvdz5, dvdz6, dvdz7, dvdz8)

    return dvdx, dvdy, dvdz
end

@inline function calcElemFBHourglassForce(xd, yd, zd,
                                  hourgam0, hourgam1,
                                  hourgam2, hourgam3,
                                  hourgam4, hourgam5,
                                  hourgam6, hourgam7,
                                  coefficient)
    i00 = 1
    i01 = 2
    i02 = 3
    i03 = 4
    @inbounds begin
    h00 = (
      hourgam0[i00] * xd[1] + hourgam1[i00] * xd[2] +
      hourgam2[i00] * xd[3] + hourgam3[i00] * xd[4] +
      hourgam4[i00] * xd[5] + hourgam5[i00] * xd[6] +
      hourgam6[i00] * xd[7] + hourgam7[i00] * xd[8]
    )

    h01 = (
      hourgam0[i01] * xd[1] + hourgam1[i01] * xd[2] +
      hourgam2[i01] * xd[3] + hourgam3[i01] * xd[4] +
      hourgam4[i01] * xd[5] + hourgam5[i01] * xd[6] +
      hourgam6[i01] * xd[7] + hourgam7[i01] * xd[8]
    )

    h02 = (
      hourgam0[i02] * xd[1] + hourgam1[i02] * xd[2] +
      hourgam2[i02] * xd[3] + hourgam3[i02] * xd[4] +
      hourgam4[i02] * xd[5] + hourgam5[i02] * xd[6] +
      hourgam6[i02] * xd[7] + hourgam7[i02] * xd[8]
    )

    h03 = (
      hourgam0[i03] * xd[1] + hourgam1[i03] * xd[2] +
      hourgam2[i03] * xd[3] + hourgam3[i03] * xd[4] +
      hourgam4[i03] * xd[5] + hourgam5[i03] * xd[6] +
      hourgam6[i03] * xd[7] + hourgam7[i03] * xd[8]
    )

    hgfx1 = coefficient *
     (hourgam0[i00] * h00 + hourgam0[i01] * h01 +
      hourgam0[i02] * h02 + hourgam0[i03] * h03)

    hgfx2 = coefficient *
     (hourgam1[i00] * h00 + hourgam1[i01] * h01 +
      hourgam1[i02] * h02 + hourgam1[i03] * h03)

    hgfx3 = coefficient *
     (hourgam2[i00] * h00 + hourgam2[i01] * h01 +
      hourgam2[i02] * h02 + hourgam2[i03] * h03)

    hgfx4 = coefficient *
     (hourgam3[i00] * h00 + hourgam3[i01] * h01 +
      hourgam3[i02] * h02 + hourgam3[i03] * h03)

    hgfx5 = coefficient *
     (hourgam4[i00] * h00 + hourgam4[i01] * h01 +
      hourgam4[i02] * h02 + hourgam4[i03] * h03)

    hgfx6 = coefficient *
     (hourgam5[i00] * h00 + hourgam5[i01] * h01 +
      hourgam5[i02] * h02 + hourgam5[i03] * h03)

    hgfx7 = coefficient *
     (hourgam6[i00] * h00 + hourgam6[i01] * h01 +
      hourgam6[i02] * h02 + hourgam6[i03] * h03)

    hgfx8 = coefficient *
     (hourgam7[i00] * h00 + hourgam7[i01] * h01 +
      hourgam7[i02] * h02 + hourgam7[i03] * h03)

    hgfx = SVector(hgfx1, hgfx2, hgfx3, hgfx4, hgfx5, hgfx6, hgfx7, hgfx8)

    h00 = (
      hourgam0[i00] * yd[1] + hourgam1[i00] * yd[2] +
      hourgam2[i00] * yd[3] + hourgam3[i00] * yd[4] +
      hourgam4[i00] * yd[5] + hourgam5[i00] * yd[6] +
      hourgam6[i00] * yd[7] + hourgam7[i00] * yd[8]
    )

    h01 = (
      hourgam0[i01] * yd[1] + hourgam1[i01] * yd[2] +
      hourgam2[i01] * yd[3] + hourgam3[i01] * yd[4] +
      hourgam4[i01] * yd[5] + hourgam5[i01] * yd[6] +
      hourgam6[i01] * yd[7] + hourgam7[i01] * yd[8]
    )

    h02 = (
      hourgam0[i02] * yd[1] + hourgam1[i02] * yd[2]+
      hourgam2[i02] * yd[3] + hourgam3[i02] * yd[4]+
      hourgam4[i02] * yd[5] + hourgam5[i02] * yd[6]+
      hourgam6[i02] * yd[7] + hourgam7[i02] * yd[8]
    )

    h03 = (
      hourgam0[i03] * yd[1] + hourgam1[i03] * yd[2] +
      hourgam2[i03] * yd[3] + hourgam3[i03] * yd[4] +
      hourgam4[i03] * yd[5] + hourgam5[i03] * yd[6] +
      hourgam6[i03] * yd[7] + hourgam7[i03] * yd[8]
    )


    hgfy1 = coefficient *
     (hourgam0[i00] * h00 + hourgam0[i01] * h01 +
      hourgam0[i02] * h02 + hourgam0[i03] * h03)

    hgfy2 = coefficient *
     (hourgam1[i00] * h00 + hourgam1[i01] * h01 +
      hourgam1[i02] * h02 + hourgam1[i03] * h03)

    hgfy3 = coefficient *
     (hourgam2[i00] * h00 + hourgam2[i01] * h01 +
      hourgam2[i02] * h02 + hourgam2[i03] * h03)

    hgfy4 = coefficient *
     (hourgam3[i00] * h00 + hourgam3[i01] * h01 +
      hourgam3[i02] * h02 + hourgam3[i03] * h03)

    hgfy5 = coefficient *
     (hourgam4[i00] * h00 + hourgam4[i01] * h01 +
      hourgam4[i02] * h02 + hourgam4[i03] * h03)

    hgfy6 = coefficient *
     (hourgam5[i00] * h00 + hourgam5[i01] * h01 +
      hourgam5[i02] * h02 + hourgam5[i03] * h03)

    hgfy7 = coefficient *
     (hourgam6[i00] * h00 + hourgam6[i01] * h01 +
      hourgam6[i02] * h02 + hourgam6[i03] * h03)

    hgfy8 = coefficient *
     (hourgam7[i00] * h00 + hourgam7[i01] * h01 +
      hourgam7[i02] * h02 + hourgam7[i03] * h03)

    hgfy = SVector(hgfy1, hgfy2, hgfy3, hgfy4, hgfy5, hgfy6, hgfy7, hgfy8)

    h00 = (
      hourgam0[i00] * zd[1] + hourgam1[i00] * zd[2] +
      hourgam2[i00] * zd[3] + hourgam3[i00] * zd[4] +
      hourgam4[i00] * zd[5] + hourgam5[i00] * zd[6] +
      hourgam6[i00] * zd[7] + hourgam7[i00] * zd[8]
    )

    h01 = (
      hourgam0[i01] * zd[1] + hourgam1[i01] * zd[2] +
      hourgam2[i01] * zd[3] + hourgam3[i01] * zd[4] +
      hourgam4[i01] * zd[5] + hourgam5[i01] * zd[6] +
      hourgam6[i01] * zd[7] + hourgam7[i01] * zd[8]
    )

    h02 =(
      hourgam0[i02] * zd[1] + hourgam1[i02] * zd[2]+
      hourgam2[i02] * zd[3] + hourgam3[i02] * zd[4]+
      hourgam4[i02] * zd[5] + hourgam5[i02] * zd[6]+
      hourgam6[i02] * zd[7] + hourgam7[i02] * zd[8]
    )

    h03 = (
      hourgam0[i03] * zd[1] + hourgam1[i03] * zd[2] +
      hourgam2[i03] * zd[3] + hourgam3[i03] * zd[4] +
      hourgam4[i03] * zd[5] + hourgam5[i03] * zd[6] +
      hourgam6[i03] * zd[7] + hourgam7[i03] * zd[8]
    )


    hgfz1 = coefficient *
     (hourgam0[i00] * h00 + hourgam0[i01] * h01 +
      hourgam0[i02] * h02 + hourgam0[i03] * h03)

    hgfz2 = coefficient *
     (hourgam1[i00] * h00 + hourgam1[i01] * h01 +
      hourgam1[i02] * h02 + hourgam1[i03] * h03)

    hgfz3 = coefficient *
     (hourgam2[i00] * h00 + hourgam2[i01] * h01 +
      hourgam2[i02] * h02 + hourgam2[i03] * h03)

    hgfz4 = coefficient *
     (hourgam3[i00] * h00 + hourgam3[i01] * h01 +
      hourgam3[i02] * h02 + hourgam3[i03] * h03)

    hgfz5 = coefficient *
     (hourgam4[i00] * h00 + hourgam4[i01] * h01 +
      hourgam4[i02] * h02 + hourgam4[i03] * h03)

    hgfz6 = coefficient *
     (hourgam5[i00] * h00 + hourgam5[i01] * h01 +
      hourgam5[i02] * h02 + hourgam5[i03] * h03)

    hgfz7 = coefficient *
     (hourgam6[i00] * h00 + hourgam6[i01] * h01 +
      hourgam6[i02] * h02 + hourgam6[i03] * h03)

    hgfz8 = coefficient *
     (hourgam7[i00] * h00 + hourgam7[i01] * h01 +
      hourgam7[i02] * h02 + hourgam7[i03] * h03)

    hgfz = SVector(hgfz1, hgfz2, hgfz3, hgfz4, hgfz5, hgfz6, hgfz7, hgfz8)

    end # inbounds

    return hgfx, hgfy, hgfz
end


function JTcalcElemFBHourglassForce(xd, yd, zd,
                                  hourgam0, hourgam1,
                                  hourgam2, hourgam3,
                                  hourgam4, hourgam5,
                                  hourgam6, hourgam7,
                                  coefficient, i, i3)
    i00 = 1
    i01 = 2
    i02 = 3
    i03 = 4
    i0 = i3
    i1 = i3 + 1
    i2 = i3 + 2
    i_3 = i3 + 3
    i4 = i3 + 4
    i5 = i3 + 5
    i6 = i3 + 6
    i7 = i3 + 7
    @inbounds begin
        # X component
        h00 = (
          hourgam0[i00, i] * xd[i0] + hourgam1[i00, i] * xd[i1] +
          hourgam2[i00, i] * xd[i2] + hourgam3[i00, i] * xd[i_3] +
          hourgam4[i00, i] * xd[i4] + hourgam5[i00, i] * xd[i5] +
          hourgam6[i00, i] * xd[i6] + hourgam7[i00, i] * xd[i7]
        )
    
        h01 = (
          hourgam0[i01, i] * xd[i0] + hourgam1[i01, i] * xd[i1] +
          hourgam2[i01, i] * xd[i2] + hourgam3[i01, i] * xd[i_3] +
          hourgam4[i01, i] * xd[i4] + hourgam5[i01, i] * xd[i5] +
          hourgam6[i01, i] * xd[i6] + hourgam7[i01, i] * xd[i7]
        )
    
        h02 = (
          hourgam0[i02, i] * xd[i0] + hourgam1[i02, i] * xd[i1] +
          hourgam2[i02, i] * xd[i2] + hourgam3[i02, i] * xd[i_3] +
          hourgam4[i02, i] * xd[i4] + hourgam5[i02, i] * xd[i5] +
          hourgam6[i02, i] * xd[i6] + hourgam7[i02, i] * xd[i7]
        )
    
        h03 = (
          hourgam0[i03, i] * xd[i0] + hourgam1[i03, i] * xd[i1] +
          hourgam2[i03, i] * xd[i2] + hourgam3[i03, i] * xd[i_3] +
          hourgam4[i03, i] * xd[i4] + hourgam5[i03, i] * xd[i5] +
          hourgam6[i03, i] * xd[i6] + hourgam7[i03, i] * xd[i7]
        )
    
        hgfx1 = coefficient *
         (hourgam0[i00, i] * h00 + hourgam0[i01, i] * h01 +
          hourgam0[i02, i] * h02 + hourgam0[i03, i] * h03)

        # fx_elem[i3] = coefficient *
        #  (hourgam0[i00, i] * h00 + hourgam0[i01, i] * h01 +
        #   hourgam0[i02, i] * h02 + hourgam0[i03, i] * h03)
    
        hgfx2 = coefficient *
         (hourgam1[i00, i] * h00 + hourgam1[i01, i] * h01 +
          hourgam1[i02, i] * h02 + hourgam1[i03, i] * h03)
    
        hgfx3 = coefficient *
         (hourgam2[i00, i] * h00 + hourgam2[i01, i] * h01 +
          hourgam2[i02, i] * h02 + hourgam2[i03, i] * h03)
    
        hgfx4 = coefficient *
         (hourgam3[i00, i] * h00 + hourgam3[i01, i] * h01 +
          hourgam3[i02, i] * h02 + hourgam3[i03, i] * h03)
    
        hgfx5 = coefficient *
         (hourgam4[i00, i] * h00 + hourgam4[i01, i] * h01 +
          hourgam4[i02, i] * h02 + hourgam4[i03, i] * h03)
    
        hgfx6 = coefficient *
         (hourgam5[i00, i] * h00 + hourgam5[i01, i] * h01 +
          hourgam5[i02, i] * h02 + hourgam5[i03, i] * h03)
    
        hgfx7 = coefficient *
         (hourgam6[i00, i] * h00 + hourgam6[i01, i] * h01 +
          hourgam6[i02, i] * h02 + hourgam6[i03, i] * h03)
    
        hgfx8 = coefficient *
         (hourgam7[i00, i] * h00 + hourgam7[i01, i] * h01 +
          hourgam7[i02, i] * h02 + hourgam7[i03, i] * h03)
    
        # hgfx = SVector(hgfx1, hgfx2, hgfx3, hgfx4, hgfx5, hgfx6, hgfx7, hgfx8)
    
        # Y component
        h00 = (
          hourgam0[i00, i] * yd[i0] + hourgam1[i00, i] * yd[i1] +
          hourgam2[i00, i] * yd[i2] + hourgam3[i00, i] * yd[i_3] +
          hourgam4[i00, i] * yd[i4] + hourgam5[i00, i] * yd[i5] +
          hourgam6[i00, i] * yd[i6] + hourgam7[i00, i] * yd[i7]
        )
    
        h01 = (
          hourgam0[i01, i] * yd[i0] + hourgam1[i01, i] * yd[i1] +
          hourgam2[i01, i] * yd[i2] + hourgam3[i01, i] * yd[i_3] +
          hourgam4[i01, i] * yd[i4] + hourgam5[i01, i] * yd[i5] +
          hourgam6[i01, i] * yd[i6] + hourgam7[i01, i] * yd[i7]
        )
    
        h02 = (
          hourgam0[i02, i] * yd[i0] + hourgam1[i02, i] * yd[i1] +
          hourgam2[i02, i] * yd[i2] + hourgam3[i02, i] * yd[i_3] +
          hourgam4[i02, i] * yd[i4] + hourgam5[i02, i] * yd[i5] +
          hourgam6[i02, i] * yd[i6] + hourgam7[i02, i] * yd[i7]
        )
    
        h03 = (
          hourgam0[i03, i] * yd[i0] + hourgam1[i03, i] * yd[i1] +
          hourgam2[i03, i] * yd[i2] + hourgam3[i03, i] * yd[i_3] +
          hourgam4[i03, i] * yd[i4] + hourgam5[i03, i] * yd[i5] +
          hourgam6[i03, i] * yd[i6] + hourgam7[i03, i] * yd[i7]
        )
    
        hgfy1 = coefficient *
         (hourgam0[i00, i] * h00 + hourgam0[i01, i] * h01 +
          hourgam0[i02, i] * h02 + hourgam0[i03, i] * h03)
    
        hgfy2 = coefficient *
         (hourgam1[i00, i] * h00 + hourgam1[i01, i] * h01 +
          hourgam1[i02, i] * h02 + hourgam1[i03, i] * h03)
    
        hgfy3 = coefficient *
         (hourgam2[i00, i] * h00 + hourgam2[i01, i] * h01 +
          hourgam2[i02, i] * h02 + hourgam2[i03, i] * h03)
    
        hgfy4 = coefficient *
         (hourgam3[i00, i] * h00 + hourgam3[i01, i] * h01 +
          hourgam3[i02, i] * h02 + hourgam3[i03, i] * h03)
    
        hgfy5 = coefficient *
         (hourgam4[i00, i] * h00 + hourgam4[i01, i] * h01 +
          hourgam4[i02, i] * h02 + hourgam4[i03, i] * h03)
    
        hgfy6 = coefficient *
         (hourgam5[i00, i] * h00 + hourgam5[i01, i] * h01 +
          hourgam5[i02, i] * h02 + hourgam5[i03, i] * h03)
    
        hgfy7 = coefficient *
         (hourgam6[i00, i] * h00 + hourgam6[i01, i] * h01 +
          hourgam6[i02, i] * h02 + hourgam6[i03, i] * h03)
    
        hgfy8 = coefficient *
         (hourgam7[i00, i] * h00 + hourgam7[i01, i] * h01 +
          hourgam7[i02, i] * h02 + hourgam7[i03, i] * h03)
    
        # hgfy = SVector(hgfy1, hgfy2, hgfy3, hgfy4, hgfy5, hgfy6, hgfy7, hgfy8)
    
        # Z component
        h00 = (
          hourgam0[i00, i] * zd[i0] + hourgam1[i00, i] * zd[i1] +
          hourgam2[i00, i] * zd[i2] + hourgam3[i00, i] * zd[i_3] +
          hourgam4[i00, i] * zd[i4] + hourgam5[i00, i] * zd[i5] +
          hourgam6[i00, i] * zd[i6] + hourgam7[i00, i] * zd[i7]
        )
    
        h01 = (
          hourgam0[i01, i] * zd[i0] + hourgam1[i01, i] * zd[i1] +
          hourgam2[i01, i] * zd[i2] + hourgam3[i01, i] * zd[i_3] +
          hourgam4[i01, i] * zd[i4] + hourgam5[i01, i] * zd[i5] +
          hourgam6[i01, i] * zd[i6] + hourgam7[i01, i] * zd[i7]
        )
    
        h02 = (
          hourgam0[i02, i] * zd[i0] + hourgam1[i02, i] * zd[i1] +
          hourgam2[i02, i] * zd[i2] + hourgam3[i02, i] * zd[i_3] +
          hourgam4[i02, i] * zd[i4] + hourgam5[i02, i] * zd[i5] +
          hourgam6[i02, i] * zd[i6] + hourgam7[i02, i] * zd[i7]
        )
    
        h03 = (
          hourgam0[i03, i] * zd[i0] + hourgam1[i03, i] * zd[i1] +
          hourgam2[i03, i] * zd[i2] + hourgam3[i03, i] * zd[i_3] +
          hourgam4[i03, i] * zd[i4] + hourgam5[i03, i] * zd[i5] +
          hourgam6[i03, i] * zd[i6] + hourgam7[i03, i] * zd[i7]
        )
    
        hgfz1 = coefficient *
         (hourgam0[i00, i] * h00 + hourgam0[i01, i] * h01 +
          hourgam0[i02, i] * h02 + hourgam0[i03, i] * h03)
    
        hgfz2 = coefficient *
         (hourgam1[i00, i] * h00 + hourgam1[i01, i] * h01 +
          hourgam1[i02, i] * h02 + hourgam1[i03, i] * h03)
    
        hgfz3 = coefficient *
         (hourgam2[i00, i] * h00 + hourgam2[i01, i] * h01 +
          hourgam2[i02, i] * h02 + hourgam2[i03, i] * h03)
    
        hgfz4 = coefficient *
         (hourgam3[i00, i] * h00 + hourgam3[i01, i] * h01 +
          hourgam3[i02, i] * h02 + hourgam3[i03, i] * h03)
    
        hgfz5 = coefficient *
         (hourgam4[i00, i] * h00 + hourgam4[i01, i] * h01 +
          hourgam4[i02, i] * h02 + hourgam4[i03, i] * h03)
    
        hgfz6 = coefficient *
         (hourgam5[i00, i] * h00 + hourgam5[i01, i] * h01 +
          hourgam5[i02, i] * h02 + hourgam5[i03, i] * h03)
    
        hgfz7 = coefficient *
         (hourgam6[i00, i] * h00 + hourgam6[i01, i] * h01 +
          hourgam6[i02, i] * h02 + hourgam6[i03, i] * h03)
    
        hgfz8 = coefficient *
         (hourgam7[i00, i] * h00 + hourgam7[i01, i] * h01 +
          hourgam7[i02, i] * h02 + hourgam7[i03, i] * h03)

        # hgfz = SVector(hgfz1, hgfz2, hgfz3, hgfz4, hgfz5, hgfz6, hgfz7, hgfz8) # Change these parts perhaps

    end # inbounds

    # return hgfx, hgfy, hgfz # May have to change this return
    # fx_elem[i3] = 0.0
    # fx_elem[i3] = Float64(hgfx1)
    return hgfx1, hgfx2, hgfx3, hgfx4, hgfx5, hgfx6, hgfx7, hgfx8, hgfy1, hgfy2, hgfy3, hgfy4, hgfy5, hgfy6, hgfy7, hgfy8, hgfz1, hgfz2, hgfz3, hgfz4, hgfz5, hgfz6, hgfz7, hgfz8
    # return hgfx1, hgfx2
end


function calcFBHourglassForceForElems_A(i, determ, x8n, y8n, z8n, hourgam0, hourgam1, hourgam2, hourgam3, hourgam4, hourgam5, hourgam6, hourgam7, dvdx, dvdy, dvdz, fx_elem, fy_elem, fz_elem, xd1, yd1, zd1, gamma, ss, elemMass, nodelist, xd, yd, zd, hourg)

    @inbounds begin

        if(i > 884736)
            # println("i: ", i)
            # looperVariable[16] == 25
            return
        end
        # Section 01


        i3 = (8 * (i - 1)) + 1
        volinv = 1.0/determ[i]

        for i1 in 1:4
            
            hourmodx = 
                x8n[i3]     * gamma[1, i1] + x8n[i3 + 1] * gamma[2, i1] + 
                x8n[i3 + 2] * gamma[3, i1] + x8n[i3 + 3] * gamma[4, i1] +
                x8n[i3 + 4] * gamma[5, i1] + x8n[i3 + 5] * gamma[6, i1] +
                x8n[i3 + 6] * gamma[7, i1] + x8n[i3 + 7] * gamma[8, i1]

            if i1 == 1
                # looperVariable[16] = hourmodx
            end

            hourmody =
                y8n[i3]     * gamma[1, i1] + y8n[i3 + 1] * gamma[2, i1] +
                y8n[i3 + 2] * gamma[3, i1] + y8n[i3 + 3] * gamma[4, i1] +
                y8n[i3 + 4] * gamma[5, i1] + y8n[i3 + 5] * gamma[6, i1] +
                y8n[i3 + 6] * gamma[7, i1] + y8n[i3 + 7] * gamma[8, i1]

            hourmodz =
                z8n[i3]     * gamma[1, i1] + z8n[i3 + 1] * gamma[2, i1] +
                z8n[i3 + 2] * gamma[3, i1] + z8n[i3 + 3] * gamma[4, i1] +
                z8n[i3 + 4] * gamma[5, i1] + z8n[i3 + 5] * gamma[6, i1] +
                z8n[i3 + 6] * gamma[7, i1] + z8n[i3 + 7] * gamma[8, i1]

            hourgam0[i1, i] = gamma[1, i1] - volinv * (dvdx[i3]     * hourmodx + 
                                dvdy[i3] * hourmody + dvdz[i3] * hourmodz)

            hourgam1[i1, i] = gamma[2, i1] - volinv * (dvdx[i3 + 1] * hourmodx + 
                                dvdy[i3 + 1] * hourmody + dvdz[i3 + 1] * hourmodz)

            hourgam2[i1, i] = gamma[3, i1] - volinv * (dvdx[i3 + 2] * hourmodx + 
                                dvdy[i3 + 2] * hourmody + dvdz[i3 + 2] * hourmodz)

            hourgam3[i1, i] = gamma[4, i1] - volinv * (dvdx[i3 + 3] * hourmodx + 
                                dvdy[i3 + 3] * hourmody + dvdz[i3 + 3] * hourmodz)

            hourgam4[i1, i] = gamma[5, i1] - volinv * (dvdx[i3 + 4] * hourmodx + 
                                dvdy[i3 + 4] * hourmody + dvdz[i3 + 4] * hourmodz)

            hourgam5[i1, i] = gamma[6, i1] - volinv * (dvdx[i3 + 5] * hourmodx + 
                                dvdy[i3 + 5] * hourmody + dvdz[i3 + 5] * hourmodz)

            hourgam6[i1, i] = gamma[7, i1] - volinv * (dvdx[i3 + 6] * hourmodx + 
                                dvdy[i3 + 6] * hourmody + dvdz[i3 + 6] * hourmodz)

            hourgam7[i1, i] = gamma[8, i1] - volinv * (dvdx[i3 + 7] * hourmodx + 
                                dvdy[i3 + 7] * hourmody + dvdz[i3 + 7] * hourmodz)

        end
        # # Section 02
        # #   compute forces
        # #   store forces into h arrays (force arrays)

        ss1 = ss[i]
        mass1 = elemMass[i]
        volume13 = cbrt(determ[i])

        n0si2 = nodelist[((i - 1) * 8) + 1]
        n1si2 = nodelist[((i - 1) * 8) + 2]
        n2si2 = nodelist[((i - 1) * 8) + 3]
        n3si2 = nodelist[((i - 1) * 8) + 4]
        n4si2 = nodelist[((i - 1) * 8) + 5]
        n5si2 = nodelist[((i - 1) * 8) + 6]
        n6si2 = nodelist[((i - 1) * 8) + 7]
        n7si2 = nodelist[((i - 1) * 8) + 8]

        xd1[i3] =     xd[n0si2 + 1]
        xd1[i3 + 1] = xd[n1si2 + 1]
        xd1[i3 + 2] = xd[n2si2 + 1]
        xd1[i3 + 3] = xd[n3si2 + 1]
        xd1[i3 + 4] = xd[n4si2 + 1]
        xd1[i3 + 5] = xd[n5si2 + 1]
        xd1[i3 + 6] = xd[n6si2 + 1]
        xd1[i3 + 7] = xd[n7si2 + 1]

        yd1[i3] =     yd[n0si2 + 1]
        yd1[i3 + 1] = yd[n1si2 + 1]
        yd1[i3 + 2] = yd[n2si2 + 1]
        yd1[i3 + 3] = yd[n3si2 + 1]
        yd1[i3 + 4] = yd[n4si2 + 1]
        yd1[i3 + 5] = yd[n5si2 + 1]
        yd1[i3 + 6] = yd[n6si2 + 1]
        yd1[i3 + 7] = yd[n7si2 + 1]

        zd1[i3] =     zd[n0si2 + 1]
        zd1[i3 + 1] = zd[n1si2 + 1]
        zd1[i3 + 2] = zd[n2si2 + 1]
        zd1[i3 + 3] = zd[n3si2 + 1]
        zd1[i3 + 4] = zd[n4si2 + 1]
        zd1[i3 + 5] = zd[n5si2 + 1]
        zd1[i3 + 6] = zd[n6si2 + 1]
        zd1[i3 + 7] = zd[n7si2 + 1]

        coefficient = - hourg * 0.01 * ss1 * mass1 / volume13
        if i == 1
            # looperVariable[1] = coefficient
            # looperVariable[2] = hourg
            # looperVariable[3] = ss1
            # looperVariable[4] = mass1
            # looperVariable[5] = volume13
            # looperVariable[6] = determ[i]
        end
        # # Section 03
        hgfx1, hgfx2, hgfx3, hgfx4, hgfx5, hgfx6, hgfx7, hgfx8, hgfy1, hgfy2, hgfy3, hgfy4, hgfy5, hgfy6, hgfy7, hgfy8, hgfz1, hgfz2, hgfz3, hgfz4, hgfz5, hgfz6, hgfz7, hgfz8 = JTcalcElemFBHourglassForce(xd1, yd1, zd1, hourgam0, hourgam1, hourgam2, hourgam3, hourgam4, hourgam5, hourgam6, hourgam7, coefficient, i, i3)
        # JTcalcElemFBHourglassForce(xd1, yd1, zd1, hourgam0, hourgam1, hourgam2, hourgam3, hourgam4, hourgam5, hourgam6, hourgam7, coefficient, i, i3, hgfx, hgfy, hgfz)
        # hgfx1, hgfx2 = JTcalcElemFBHourglassForce(xd1, yd1, zd1, hourgam0, hourgam1, hourgam2, hourgam3, hourgam4, hourgam5, hourgam6, hourgam7, coefficient, i, i3)

        @inbounds begin

            # if(i3 > 7077881)
            #     return
            # end
            # if(i > 884736)
            #     return
            # end
            # if i == 1
            #     looperVariable[7] = hgfx1
            # end

            # fx_elem[i3] = 0
            # hgfx1 = 0.0
            # hgfx2 = 0.0
            # hgfx3 = 0.0
            # hgfx4 = 0.0
            # hgfx5 = 0.0
            # hgfx6 = 0.0
            # hgfx7 = 0.0
            # hgfx8 = 0.0

            # fx_elem[i3] =     0.0
            fx_elem[i3] =     hgfx1
            fx_elem[i3 + 1] = hgfx2
            fx_elem[i3 + 2] = hgfx3
            fx_elem[i3 + 3] = hgfx4
            fx_elem[i3 + 4] = hgfx5
            fx_elem[i3 + 5] = hgfx6
            fx_elem[i3 + 6] = hgfx7
            fx_elem[i3 + 7] = hgfx8

            fy_elem[i3] =     hgfy1
            fy_elem[i3 + 1] = hgfy2
            fy_elem[i3 + 2] = hgfy3
            fy_elem[i3 + 3] = hgfy4
            fy_elem[i3 + 4] = hgfy5
            fy_elem[i3 + 5] = hgfy6
            fy_elem[i3 + 6] = hgfy7
            fy_elem[i3 + 7] = hgfy8

            fz_elem[i3] =     hgfz1
            fz_elem[i3 + 1] = hgfz2
            fz_elem[i3 + 2] = hgfz3
            fz_elem[i3 + 3] = hgfz4
            fz_elem[i3 + 4] = hgfz5
            fz_elem[i3 + 5] = hgfz6
            fz_elem[i3 + 6] = hgfz7
            fz_elem[i3 + 7] = hgfz8

            # fx_elem[i3 + 1] = 0.0
            # fx_elem[i3 + 2] = 0.0
            # fx_elem[i3 + 3] = 0.0
            # fx_elem[i3 + 4] = 0.0
            # fx_elem[i3 + 5] = 0.0
            # fx_elem[i3 + 6] = 0.0
            # fx_elem[i3 + 7] = 0.0



            # if i == 1
            #     for loops in 1:8
            #         hgfx[loops] = fx_elem[loops]
            #         hgfy[loops] = fy_elem[loops]
            #         # hgfz[loops] = fz_elem[loops]
            #     end
            #     looperVariable[1] = coefficient
            #     # looperVariable[2] = hourg
            #     # looperVariable[3] = ss1
            #     # looperVariable[4] = mass1
            #     # looperVariable[5] = volume13
            #     # looperVariable[6] = determ[i]

            # end

        end

    end

end


function MultiDcalcFBHourglassForceForElems(domain, determ,
                                        x8n, y8n, z8n,
                                        dvdx, dvdy, dvdz,
                                        hourg             )

    # *************************************************
    # *
    # *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
    # *               force.
    # *
    # *************************************************

    numElem = domain.numElem
    numElem8 = numElem * 8
    # println("numElem8: ", numElem8)

    fx_elem = Vector{Float64}(undef, numElem8) # Multidimensionals
    fy_elem = Vector{Float64}(undef, numElem8) 
    fz_elem = Vector{Float64}(undef, numElem8)
    # fx_elem = zeros(Float64, numElem8)



    xd1 = Vector{Float64}(undef, numElem8)
    yd1 = Vector{Float64}(undef, numElem8)
    zd1 = Vector{Float64}(undef, numElem8)

    hourgam0 = MVector{4, Float64}(undef) # Multidimensionals
    hourgam1 = MVector{4, Float64}(undef)
    hourgam2 = MVector{4, Float64}(undef)
    hourgam3 = MVector{4, Float64}(undef)
    hourgam4 = MVector{4, Float64}(undef)
    hourgam5 = MVector{4, Float64}(undef)
    hourgam6 = MVector{4, Float64}(undef)
    hourgam7 = MVector{4, Float64}(undef)


    gamma = @SMatrix [
        1.0   1.0   1.0  -1.0
        1.0  -1.0  -1.0   1.0
       -1.0  -1.0   1.0  -1.0
       -1.0   1.0  -1.0   1.0
       -1.0  -1.0   1.0   1.0
       -1.0   1.0  -1.0  -1.0
        1.0   1.0   1.0   1.0
        1.0  -1.0  -1.0  -1.0
    ]
    

#     gamma = Array{Float64}(undef, 8, 4)

# # Assign values element by element
#     gamma[1, 1] = 1.0
#     gamma[1, 2] = 1.0
#     gamma[1, 3] = 1.0
#     gamma[1, 4] = -1.0

#     gamma[2, 1] = 1.0
#     gamma[2, 2] = -1.0
#     gamma[2, 3] = -1.0
#     gamma[2, 4] = 1.0

#     gamma[3, 1] = -1.0
#     gamma[3, 2] = -1.0
#     gamma[3, 3] = 1.0
#     gamma[3, 4] = -1.0

#     gamma[4, 1] = -1.0
#     gamma[4, 2] = 1.0
#     gamma[4, 3] = -1.0
#     gamma[4, 4] = 1.0

#     gamma[5, 1] = -1.0
#     gamma[5, 2] = -1.0
#     gamma[5, 3] = 1.0
#     gamma[5, 4] = 1.0

#     gamma[6, 1] = -1.0
#     gamma[6, 2] = 1.0
#     gamma[6, 3] = -1.0
#     gamma[6, 4] = -1.0

#     gamma[7, 1] = 1.0
#     gamma[7, 2] = 1.0
#     gamma[7, 3] = 1.0
#     gamma[7, 4] = 1.0

#     gamma[8, 1] = 1.0
#     gamma[8, 2] = -1.0
#     gamma[8, 3] = -1.0
#     gamma[8, 4] = -1.0


    # fx_elem_domainVariables = JACC.Array(zeros(Float64, numElem8, 1000))
    # fy_elem_domainVariables = JACC.Array(zeros(Float64, numElem8, 1000))
    # fz_elem_domainVariables = JACC.Array(zeros(Float64, numElem8, 1000))
    hourgam0_domainVariables = JACC.Array(zeros(Float64, 4, numElem))
    hourgam1_domainVariables = JACC.Array(zeros(Float64, 4, numElem))
    hourgam2_domainVariables = JACC.Array(zeros(Float64, 4, numElem))
    hourgam3_domainVariables = JACC.Array(zeros(Float64, 4, numElem))
    hourgam4_domainVariables = JACC.Array(zeros(Float64, 4, numElem))
    hourgam5_domainVariables = JACC.Array(zeros(Float64, 4, numElem))
    hourgam6_domainVariables = JACC.Array(zeros(Float64, 4, numElem))
    hourgam7_domainVariables = JACC.Array(zeros(Float64, 4, numElem))

    hgfx = Vector{Float64}(undef, 8)
    hgfy = Vector{Float64}(undef, 8)
    hgfz = Vector{Float64}(undef, 8)

    # Initialize an array of 20 indexes with all values being 0

    
    gamma_d = JACC.Array(gamma)
    x8n_d = JACC.Array(x8n)
    y8n_d = JACC.Array(y8n)
    z8n_d = JACC.Array(z8n)
    dvdx_d = JACC.Array(dvdx)
    dvdy_d = JACC.Array(dvdy)
    dvdz_d = JACC.Array(dvdz)
    determ_d = JACC.Array(determ)
    ss_d = JACC.Array(domain.ss)
    eleMass_d = JACC.Array(domain.elemMass)
    fx_elem_d = JACC.Array(fx_elem)
    fy_elem_d = JACC.Array(fy_elem)
    fz_elem_d = JACC.Array(fz_elem)
    nodelist_d = JACC.Array(domain.nodelist)
    xd_d = JACC.Array(domain.xd)
    yd_d = JACC.Array(domain.yd)
    zd_d = JACC.Array(domain.zd)
    xd1_d = JACC.Array(xd1)
    yd1_d = JACC.Array(yd1)
    zd1_d = JACC.Array(zd1)
    nodeElemCount_d = JACC.Array(domain.nodeElemCount)
    nodeElemStart_d = JACC.Array(domain.nodeElemStart)
    nodeElemCornerList_d = JACC.Array(domain.nodeElemCornerList)
    fx_d = JACC.Array(domain.fx)
    fy_d = JACC.Array(domain.fy)
    fz_d = JACC.Array(domain.fz)
    hgfx_d = JACC.Array(hgfx)
    hgfy_d = JACC.Array(hgfy)
    hgfz_d = JACC.Array(hgfz)
    randomTestWords_d = JACC.Array(zeros(Float64, 20))

    # xd1 = JACC.Array()





    # *************************************************
    # compute the hourglass modes

    # From Here 

    # for chunk_start in 1:chunk_size:numElem
    #     chunk_end = min(chunk_start + chunk_size - 1, numElem)

    #     fx_elem_domainVariables = JACC.Array(zeros(Float64, numElem8, chunk_size))
    #     fy_elem_domainVariables = JACC.Array(zeros(Float64, numElem8, chunk_size))
    #     fz_elem_domainVariables = JACC.Array(zeros(Float64, numElem8, chunk_size))
    #     f[3, 1] = 5
    #     # JACC.parallel_for(chunk_end, calcFBHourglassForceForElems_A, chunk_start, determ_d, x8n_d, y8n_d, z8n_d, hourgam0_domainVariables, hourgam1_domainVariables, hourgam2_domainVariables, hourgam3_domainVariables, hourgam4_domainVariables, hourgam5_domainVariables, hourgam6_domainVariables, hourgam7_domainVariables, dvdx_d, dvdy_d, dvdz_d)


    # end
    
    # for i in entireLoop

    #     create array of size 8

    #     call jacc.parallel (i , neededMultiArrayButNormal, pass_array_of_8)


    #         array_of_8[1] = neededMultiArrayButNormal[i]
    #         array_of_8[2] = neededMultiArrayButNormal[i+1]
    # end
    # println("Volinv: ", 1.0 / determ[1])
    # volinv= 1.0/determ[1]
    # println("NumElem: ", numElem)

    @assert length(fx_elem_d) == numElem8 "Memory allocation failed" 
    # time1 = time()
    # JACC.parallel_for(numElem, calcFBHourglassForceForElems_A, determ_d, x8n_d, y8n_d, z8n_d, hourgam0_domainVariables, hourgam1_domainVariables, hourgam2_domainVariables, hourgam3_domainVariables, hourgam4_domainVariables, hourgam5_domainVariables, hourgam6_domainVariables, hourgam7_domainVariables, dvdx_d, dvdy_d, dvdz_d, fx_elem_d, fy_elem_d, fz_elem_d, xd1_d, yd1_d, zd1_d, gamma_d, ss_d, eleMass_d, nodelist_d, xd_d, yd_d, zd_d, hourg)
    # # JACC.parallel_for(numElem, calcFBHourglassForceForElems_A, determ_d, x8n_d, y8n_d, z8n_d, hourgam0_domainVariables, hourgam1_domainVariables, hourgam2_domainVariables, hourgam3_domainVariables, hourgam4_domainVariables, hourgam5_domainVariables, hourgam6_domainVariables, hourgam7_domainVariables, dvdx_d, dvdy_d, dvdz_d, fx_elem_d, fy_elem_d, fz_elem_d, xd1_d, yd1_d, zd1_d, gamma_d, ss_d, eleMass_d, nodelist_d, xd_d, yd_d, zd_d, hourg, hgfx_d, hgfy_d, randomTestWords_d)
    # time2 = time()
    # println("Time taken for the calcFBHourglassForceForElems_A: ", time2 - time1)

    numNode = domain.numNode

    time1 = time()
    # Perform a warmup run
    JACC.parallel_for(numElem, calcFBHourglassForceForElems_A, determ_d, x8n_d, y8n_d, z8n_d, hourgam0_domainVariables, hourgam1_domainVariables, hourgam2_domainVariables, hourgam3_domainVariables, hourgam4_domainVariables, hourgam5_domainVariables, hourgam6_domainVariables, hourgam7_domainVariables, dvdx_d, dvdy_d, dvdz_d, fx_elem_d, fy_elem_d, fz_elem_d, xd1_d, yd1_d, zd1_d, gamma_d, ss_d, eleMass_d, nodelist_d, xd_d, yd_d, zd_d, hourg)
    JACC.parallel_for(numNode, Jupdate_domain, fx_d, fy_d, fz_d, nodeElemCount_d, nodeElemStart_d, nodeElemCornerList_d, fx_elem_d, fy_elem_d, fz_elem_d)
    # CUDA.synchronize()
    time2 = time()
    println("Time taken for the calcFBHourglassForceForElems_A: ", time2 - time1)


    
    # Measure the time for 10 executions and calculate the average time
    total_time = 0.0
    num_runs = 10
    # numElem = 32 * 32 * 32
    for i in 1:num_runs
        time1 = time()
        JACC.parallel_for(numElem, calcFBHourglassForceForElems_A, determ_d, x8n_d, y8n_d, z8n_d, hourgam0_domainVariables, hourgam1_domainVariables, hourgam2_domainVariables, hourgam3_domainVariables, hourgam4_domainVariables, hourgam5_domainVariables, hourgam6_domainVariables, hourgam7_domainVariables, dvdx_d, dvdy_d, dvdz_d, fx_elem_d, fy_elem_d, fz_elem_d, xd1_d, yd1_d, zd1_d, gamma_d, ss_d, eleMass_d, nodelist_d, xd_d, yd_d, zd_d, hourg)
        JACC.parallel_for(numNode, Jupdate_domain, fx_d, fy_d, fz_d, nodeElemCount_d, nodeElemStart_d, nodeElemCornerList_d, fx_elem_d, fy_elem_d, fz_elem_d)

        time2 = time()
        total_time += (time2 - time1)
    end

    average_time = total_time / num_runs
    println("Average time taken for the calcFBHourglassForceForElems_A over $num_runs runs: ", average_time)

    # hgfx = Array(hgfx_d)
    # hgfy = Array(hgfy_d)
    # hgfz = JACC.Array(hgfz_d)
    randomTestWords = Array(randomTestWords_d)
    
    # Checking methods
    # println("hgfx: ", hgfx)
    # println("hgfy: ", hgfy)

    # println("randomTestWords: ", randomTestWords)
    # println("coefficient: ", randomTestWords[1])
    # println("hourg: ", randomTestWords[2])
    # println("ss1: ", randomTestWords[3])
    # println("mass1: ", randomTestWords[4])
    # println("volume13: ", randomTestWords[5])
    # println("determ: ", randomTestWords[6])

    # println("hgfx1: ", randomTestWords[7])
    # println("hourmodx: ", randomTestWords[16])
    # xd1_d = Array(xd1_d)
    
    # for jh in 1:8
    #     println("xd1_d: ", xd1_d[jh])
    # end


    # println("hgfx: ", hgfx)
    # println("hgfy: ", hgfy)
    # println("hgfz: ", hgfz)

    # @inbounds for i2 in 1:numElem

    #     i3=8*(i2-1)+1
    #     volinv= 1.0/determ[i2]

    #     for i1 in 1:4

    #         hourmodx =
    #             x8n[i3]   * gamma[1,i1] + x8n[i3+1] * gamma[2,i1] +
    #             x8n[i3+2] * gamma[3,i1] + x8n[i3+3] * gamma[4,i1] +
    #             x8n[i3+4] * gamma[5,i1] + x8n[i3+5] * gamma[6,i1] +
    #             x8n[i3+6] * gamma[7,i1] + x8n[i3+7] * gamma[8,i1]

    #         hourmody =
    #             y8n[i3]   * gamma[1,i1] + y8n[i3+1] * gamma[2,i1] +
    #             y8n[i3+2] * gamma[3,i1] + y8n[i3+3] * gamma[4,i1] +
    #             y8n[i3+4] * gamma[5,i1] + y8n[i3+5] * gamma[6,i1] +
    #             y8n[i3+6] * gamma[7,i1] + y8n[i3+7] * gamma[8,i1]

    #         hourmodz =
    #             z8n[i3]   * gamma[1,i1] + z8n[i3+1] * gamma[2,i1] +
    #             z8n[i3+2] * gamma[3,i1] + z8n[i3+3] * gamma[4,i1] +
    #             z8n[i3+4] * gamma[5,i1] + z8n[i3+5] * gamma[6,i1] +
    #             z8n[i3+6] * gamma[7,i1] + z8n[i3+7] * gamma[8,i1]

    #         hourgam0[i1] = gamma[1,i1] -  volinv*(dvdx[i3  ] * hourmodx +
    #                         dvdy[i3  ] * hourmody + dvdz[i3  ] * hourmodz )

    #         hourgam1[i1] = gamma[2,i1] -  volinv*(dvdx[i3+1] * hourmodx +
    #                         dvdy[i3+1] * hourmody + dvdz[i3+1] * hourmodz )

    #         hourgam2[i1] = gamma[3,i1] -  volinv*(dvdx[i3+2] * hourmodx +
    #                         dvdy[i3+2] * hourmody + dvdz[i3+2] * hourmodz )

    #         hourgam3[i1] = gamma[4,i1] -  volinv*(dvdx[i3+3] * hourmodx +
    #                         dvdy[i3+3] * hourmody + dvdz[i3+3] * hourmodz )

    #         hourgam4[i1] = gamma[5,i1] -  volinv*(dvdx[i3+4] * hourmodx +
    #                         dvdy[i3+4] * hourmody + dvdz[i3+4] * hourmodz )

    #         hourgam5[i1] = gamma[6,i1] -  volinv*(dvdx[i3+5] * hourmodx +
    #                         dvdy[i3+5] * hourmody + dvdz[i3+5] * hourmodz )

    #         hourgam6[i1] = gamma[7,i1] -  volinv*(dvdx[i3+6] * hourmodx +
    #                         dvdy[i3+6] * hourmody + dvdz[i3+6] * hourmodz )

    #         hourgam7[i1] = gamma[8,i1] -  volinv*(dvdx[i3+7] * hourmodx +
    #                         dvdy[i3+7] * hourmody + dvdz[i3+7] * hourmodz )
    #     end


    #     #   compute forces
    #     #   store forces into h arrays (force arrays)

    #     # Second Start

    #     ss1 = domain.ss[i2]
    #     mass1 = domain.elemMass[i2]
    #     volume13=cbrt(determ[i2])

    #     n0si2 = domain.nodelist[(i2-1)*8+1]
    #     n1si2 = domain.nodelist[(i2-1)*8+2]
    #     n2si2 = domain.nodelist[(i2-1)*8+3]
    #     n3si2 = domain.nodelist[(i2-1)*8+4]
    #     n4si2 = domain.nodelist[(i2-1)*8+5]
    #     n5si2 = domain.nodelist[(i2-1)*8+6]
    #     n6si2 = domain.nodelist[(i2-1)*8+7]
    #     n7si2 = domain.nodelist[(i2-1)*8+8]

    #     xd1 = SVector(
    #         domain.xd[n0si2+1], 
    #         domain.xd[n1si2+1],
    #         domain.xd[n2si2+1],
    #         domain.xd[n3si2+1],
    #         domain.xd[n4si2+1],
    #         domain.xd[n5si2+1],
    #         domain.xd[n6si2+1],
    #         domain.xd[n7si2+1],
    #     ) # Create these outside of the JACC Loop
 
    #     yd1 = SVector(
    #         domain.yd[n0si2+1],
    #         domain.yd[n1si2+1],
    #         domain.yd[n2si2+1],
    #         domain.yd[n3si2+1],
    #         domain.yd[n4si2+1],
    #         domain.yd[n5si2+1],
    #         domain.yd[n6si2+1],
    #         domain.yd[n7si2+1],
    #     )

    #     zd1 = SVector(
    #         domain.zd[n0si2+1],
    #         domain.zd[n1si2+1],
    #         domain.zd[n2si2+1],
    #         domain.zd[n3si2+1],
    #         domain.zd[n4si2+1],
    #         domain.zd[n5si2+1],
    #         domain.zd[n6si2+1],
    #         domain.zd[n7si2+1],
    #     )

    #     coefficient = - hourg * 0.01 * ss1 * mass1 / volume13

    #     # Next Pause

    #     # Second Start

    #     hgfx, hgfy, hgfz = calcElemFBHourglassForce(xd1,yd1,zd1,
    #                              hourgam0,hourgam1,hourgam2,hourgam3,
    #                              hourgam4,hourgam5,hourgam6,hourgam7,
    #                              coefficient)

    #     @inbounds begin
    #         fx_elem[i3] = hgfx[1] # Multidimensionals
    #         fx_elem[i3+1] = hgfx[2]
    #         fx_elem[i3+2] = hgfx[3]
    #         fx_elem[i3+3] = hgfx[4]
    #         fx_elem[i3+4] = hgfx[5]
    #         fx_elem[i3+5] = hgfx[6]
    #         fx_elem[i3+6] = hgfx[7]
    #         fx_elem[i3+7] = hgfx[8]

    #         fy_elem[i3] = hgfy[1]
    #         fy_elem[i3+1] = hgfy[2]
    #         fy_elem[i3+2] = hgfy[3]
    #         fy_elem[i3+3] = hgfy[4]
    #         fy_elem[i3+4] = hgfy[5]
    #         fy_elem[i3+5] = hgfy[6]
    #         fy_elem[i3+6] = hgfy[7]
    #         fy_elem[i3+7] = hgfy[8]

    #         fz_elem[i3] = hgfz[1]
    #         fz_elem[i3+1] = hgfz[2]
    #         fz_elem[i3+2] = hgfz[3]
    #         fz_elem[i3+3] = hgfz[4]
    #         fz_elem[i3+4] = hgfz[5]
    #         fz_elem[i3+5] = hgfz[6]
    #         fz_elem[i3+6] = hgfz[7]
    #         fz_elem[i3+7] = hgfz[8]
    #     end

    #     # #if 0
    #     #     domain%m_fx(n0si2) = domain%m_fx(n0si2) + hgfx(0)
    #     #     domain%m_fy(n0si2) = domain%m_fy(n0si2) + hgfy(0)
    #     #     domain%m_fz(n0si2) = domain%m_fz(n0si2) + hgfz(0)

    #     #     domain%m_fx(n1si2) = domain%m_fx(n1si2) + hgfx(1)
    #     #     domain%m_fy(n1si2) = domain%m_fy(n1si2) + hgfy(1)
    #     #     domain%m_fz(n1si2) = domain%m_fz(n1si2) + hgfz(1)

    #     #     domain%m_fx(n2si2) = domain%m_fx(n2si2) + hgfx(2)
    #     #     domain%m_fy(n2si2) = domain%m_fy(n2si2) + hgfy(2)
    #     #     domain%m_fz(n2si2) = domain%m_fz(n2si2) + hgfz(2)

    #     #     domain%m_fx(n3si2) = domain%m_fx(n3si2) + hgfx(3)
    #     #     domain%m_fy(n3si2) = domain%m_fy(n3si2) + hgfy(3)
    #     #     domain%m_fz(n3si2) = domain%m_fz(n3si2) + hgfz(3)

    #     #     domain%m_fx(n4si2) = domain%m_fx(n4si2) + hgfx(4)
    #     #     domain%m_fy(n4si2) = domain%m_fy(n4si2) + hgfy(4)
    #     #     domain%m_fz(n4si2) = domain%m_fz(n4si2) + hgfz(4)

    #     #     domain%m_fx(n5si2) = domain%m_fx(n5si2) + hgfx(5)
    #     #     domain%m_fy(n5si2) = domain%m_fy(n5si2) + hgfy(5)
    #     #     domain%m_fz(n5si2) = domain%m_fz(n5si2) + hgfz(5)

    #     #     domain%m_fx(n6si2) = domain%m_fx(n6si2) + hgfx(6)
    #     #     domain%m_fy(n6si2) = domain%m_fy(n6si2) + hgfy(6)
    #     #     domain%m_fz(n6si2) = domain%m_fz(n6si2) + hgfz(6)

    #     #     domain%m_fx(n7si2) = domain%m_fx(n7si2) + hgfx(7)
    #     #     domain%m_fy(n7si2) = domain%m_fy(n7si2) + hgfy(7)
    #     #     domain%m_fz(n7si2) = domain%m_fz(n7si2) + hgfz(7)
    #     # #endif

    #     # End Loop here
    # end

    #Here
    
   

    



    # xd = Array(xd_d)
    # copyto!(fx_elem, fx_elem_d)
    # copyto!(fy_elem, fy_elem_d)
    # copyto!(fz_elem, fz_elem_d)
    # copyto!(determ, determ_d)

    # for jh in 1:8
    #     println("fx_elem[", jh, "]: ", fx_elem[jh])
    #     println("fy_elem[", jh, "]: ", fy_elem[jh])
    #     println("fz_elem[", jh, "]: ", fz_elem[jh])
    # end

    # @inbounds for gnode in 1:numNode
    #     count = domain.nodeElemCount[gnode]
    #     start = domain.nodeElemStart[gnode]
    #     fx = 0.0
    #     fy = 0.0
    #     fz = 0.0
    #     @simd for i in 1:count
    #         elem = domain.nodeElemCornerList[start+i]
    #         fx += fx_elem[elem]
    #         fy += fy_elem[elem]
    #         fz += fz_elem[elem]
    #     end
    #     domain.fx[gnode] += fx
    #     domain.fy[gnode] += fy 
    #     domain.fz[gnode] += fz
    # end

    # resultsX = Vector{Float64}(undef, 1)
    # resultsY = Vector{Float64}(undef, 1)
    # resultsZ = Vector{Float64}(undef, 1)
    # results = Vector{Float64}(undef, 3)
    # results = zeroes(Float64, 3)
    # resultsX_d = JACC.Array(resultsX)
    # resultsY_d = JACC.Array(resultsY)
    # resultsZ_d = JACC.Array(resultsZ)

    numNode = domain.numNode



    # zeroAdder = JACC.Array(zeros(Float64, numElem8))
    # println("numNode: ", numNode)

    # fx_elem = Array(fx_elem_d)
    # fy_elem = Array(fy_elem_d)
    # fz_elem = Array(fz_elem_d)
    # @inbounds for gnode in 1:numNode
    #     count = domain.nodeElemCount[gnode]
    #     start = domain.nodeElemStart[gnode]
    #     fx = 0.0
    #     fy = 0.0
    #     fz = 0.0
    #     @simd for i in 1:count
    #         elem = domain.nodeElemCornerList[start+i]
    #         fx += fx_elem[elem]
    #         fy += fy_elem[elem]
    #         fz += fz_elem[elem]
    #     end
    #     # results_d = JACC.parallel_reduce(count, reduce_func, results_d, start, nodeElemCornerList_d, fx_elem_d, fy_elem_d, fz_elem_d)
        # resultsX_d = JACC.parallel_reduce(count, reduce_func, zeroAdder, start, nodeElemCornerList_d, fx_elem_d)
        # resultsY_d = JACC.parallel_reduce(count, reduce_func, zeroAdder, start, nodeElemCornerList_d, fy_elem_d)
        # resultsZ_d = JACC.parallel_reduce(count, reduce_func, zeroAdder, start, nodeElemCornerList_d, fz_elem_d)
    #     # resultsX = Array(resultsX_d)
    #     # resultsY = Array(resultsY_d)
    #     # resultsZ = Array(resultsZ_d)
    #     domain.fx[gnode] += fx
    #     domain.fy[gnode] += fy 
    #     domain.fz[gnode] += fz
    #     # domain.fx[gnode] += resultsX[1]
    #     # domain.fy[gnode] += resultsY[1]
    #     # domain.fz[gnode] += resultsZ[1]
    # end
    # JACC.parallel_for(numNode, Jminiupdate_domain, fx_d, fy_d, fz_d, nodeElemCount_d, nodeElemStart_d, nodeElemCornerList_d, fx_elem_d, fy_elem_d, fz_elem_d, resultsX_d, resultsY_d, resultsZ_d)   
    # JACC.parallel_for(numNode, Jupdate_domain, fx_d, fy_d, fz_d, nodeElemCount_d, nodeElemStart_d, nodeElemCornerList_d, fx_elem_d, fy_elem_d, fz_elem_d)
    
    domain.fx = Array(fx_d)
    domain.fy = Array(fy_d)
    domain.fz = Array(fz_d)

    # copyto!(domain.fx, fx_d)
    # copyto!(domain.fy, fy_d)
    # copyto!(domain.fz, fz_d)

end


@inline function calcFBHourglassForceForElems(domain, determ,
                                        x8n, y8n, z8n,
                                        dvdx, dvdy, dvdz,
                                        hourg             )

    # *************************************************
    # *
    # *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
    # *               force.
    # *
    # *************************************************

    numElem = domain.numElem
    numElem8 = numElem * 8

    fx_elem = Vector{Float64}(undef, numElem8)
    fy_elem = Vector{Float64}(undef, numElem8)
    fz_elem = Vector{Float64}(undef, numElem8)

    hourgam0 = MVector{4, Float64}(undef)
    hourgam1 = MVector{4, Float64}(undef)
    hourgam2 = MVector{4, Float64}(undef)
    hourgam3 = MVector{4, Float64}(undef)
    hourgam4 = MVector{4, Float64}(undef)
    hourgam5 = MVector{4, Float64}(undef)
    hourgam6 = MVector{4, Float64}(undef)
    hourgam7 = MVector{4, Float64}(undef)

    gamma = @SMatrix [
        1.0   1.0   1.0  -1.0
        1.0  -1.0  -1.0   1.0
       -1.0  -1.0   1.0  -1.0
       -1.0   1.0  -1.0   1.0
       -1.0  -1.0   1.0   1.0
       -1.0   1.0  -1.0  -1.0
        1.0   1.0   1.0   1.0
        1.0  -1.0  -1.0  -1.0
    ]

    # *************************************************
    # compute the hourglass modes


    @inbounds for i2 in 1:numElem

        i3=8*(i2-1)+1
        volinv= 1.0/determ[i2]

        for i1 in 1:4

            hourmodx =
                x8n[i3]   * gamma[1,i1] + x8n[i3+1] * gamma[2,i1] +
                x8n[i3+2] * gamma[3,i1] + x8n[i3+3] * gamma[4,i1] +
                x8n[i3+4] * gamma[5,i1] + x8n[i3+5] * gamma[6,i1] +
                x8n[i3+6] * gamma[7,i1] + x8n[i3+7] * gamma[8,i1]

            hourmody =
                y8n[i3]   * gamma[1,i1] + y8n[i3+1] * gamma[2,i1] +
                y8n[i3+2] * gamma[3,i1] + y8n[i3+3] * gamma[4,i1] +
                y8n[i3+4] * gamma[5,i1] + y8n[i3+5] * gamma[6,i1] +
                y8n[i3+6] * gamma[7,i1] + y8n[i3+7] * gamma[8,i1]

            hourmodz =
                z8n[i3]   * gamma[1,i1] + z8n[i3+1] * gamma[2,i1] +
                z8n[i3+2] * gamma[3,i1] + z8n[i3+3] * gamma[4,i1] +
                z8n[i3+4] * gamma[5,i1] + z8n[i3+5] * gamma[6,i1] +
                z8n[i3+6] * gamma[7,i1] + z8n[i3+7] * gamma[8,i1]

            hourgam0[i1] = gamma[1,i1] -  volinv*(dvdx[i3  ] * hourmodx +
                            dvdy[i3  ] * hourmody + dvdz[i3  ] * hourmodz )

            hourgam1[i1] = gamma[2,i1] -  volinv*(dvdx[i3+1] * hourmodx +
                            dvdy[i3+1] * hourmody + dvdz[i3+1] * hourmodz )

            hourgam2[i1] = gamma[3,i1] -  volinv*(dvdx[i3+2] * hourmodx +
                            dvdy[i3+2] * hourmody + dvdz[i3+2] * hourmodz )

            hourgam3[i1] = gamma[4,i1] -  volinv*(dvdx[i3+3] * hourmodx +
                            dvdy[i3+3] * hourmody + dvdz[i3+3] * hourmodz )

            hourgam4[i1] = gamma[5,i1] -  volinv*(dvdx[i3+4] * hourmodx +
                            dvdy[i3+4] * hourmody + dvdz[i3+4] * hourmodz )

            hourgam5[i1] = gamma[6,i1] -  volinv*(dvdx[i3+5] * hourmodx +
                            dvdy[i3+5] * hourmody + dvdz[i3+5] * hourmodz )

            hourgam6[i1] = gamma[7,i1] -  volinv*(dvdx[i3+6] * hourmodx +
                            dvdy[i3+6] * hourmody + dvdz[i3+6] * hourmodz )

            hourgam7[i1] = gamma[8,i1] -  volinv*(dvdx[i3+7] * hourmodx +
                            dvdy[i3+7] * hourmody + dvdz[i3+7] * hourmodz )
        end


        #   compute forces
        #   store forces into h arrays (force arrays)

        ss1 = domain.ss[i2]
        mass1 = domain.elemMass[i2]
        volume13=cbrt(determ[i2])

        n0si2 = domain.nodelist[(i2-1)*8+1]
        n1si2 = domain.nodelist[(i2-1)*8+2]
        n2si2 = domain.nodelist[(i2-1)*8+3]
        n3si2 = domain.nodelist[(i2-1)*8+4]
        n4si2 = domain.nodelist[(i2-1)*8+5]
        n5si2 = domain.nodelist[(i2-1)*8+6]
        n6si2 = domain.nodelist[(i2-1)*8+7]
        n7si2 = domain.nodelist[(i2-1)*8+8]

        xd1 = SVector(
            domain.xd[n0si2+1],
            domain.xd[n1si2+1],
            domain.xd[n2si2+1],
            domain.xd[n3si2+1],
            domain.xd[n4si2+1],
            domain.xd[n5si2+1],
            domain.xd[n6si2+1],
            domain.xd[n7si2+1],
        )

        yd1 = SVector(
            domain.yd[n0si2+1],
            domain.yd[n1si2+1],
            domain.yd[n2si2+1],
            domain.yd[n3si2+1],
            domain.yd[n4si2+1],
            domain.yd[n5si2+1],
            domain.yd[n6si2+1],
            domain.yd[n7si2+1],
        )

        zd1 = SVector(
            domain.zd[n0si2+1],
            domain.zd[n1si2+1],
            domain.zd[n2si2+1],
            domain.zd[n3si2+1],
            domain.zd[n4si2+1],
            domain.zd[n5si2+1],
            domain.zd[n6si2+1],
            domain.zd[n7si2+1],
        )

        coefficient = - hourg * 0.01 * ss1 * mass1 / volume13

        hgfx, hgfy, hgfz = calcElemFBHourglassForce(xd1,yd1,zd1,
                                 hourgam0,hourgam1,hourgam2,hourgam3,
                                 hourgam4,hourgam5,hourgam6,hourgam7,
                                 coefficient)

        @inbounds begin
            fx_elem[i3] = hgfx[1]
            fx_elem[i3+1] = hgfx[2]
            fx_elem[i3+2] = hgfx[3]
            fx_elem[i3+3] = hgfx[4]
            fx_elem[i3+4] = hgfx[5]
            fx_elem[i3+5] = hgfx[6]
            fx_elem[i3+6] = hgfx[7]
            fx_elem[i3+7] = hgfx[8]

            fy_elem[i3] = hgfy[1]
            fy_elem[i3+1] = hgfy[2]
            fy_elem[i3+2] = hgfy[3]
            fy_elem[i3+3] = hgfy[4]
            fy_elem[i3+4] = hgfy[5]
            fy_elem[i3+5] = hgfy[6]
            fy_elem[i3+6] = hgfy[7]
            fy_elem[i3+7] = hgfy[8]

            fz_elem[i3] = hgfz[1]
            fz_elem[i3+1] = hgfz[2]
            fz_elem[i3+2] = hgfz[3]
            fz_elem[i3+3] = hgfz[4]
            fz_elem[i3+4] = hgfz[5]
            fz_elem[i3+5] = hgfz[6]
            fz_elem[i3+6] = hgfz[7]
            fz_elem[i3+7] = hgfz[8]
        end

    # #if 0
    #     domain%m_fx(n0si2) = domain%m_fx(n0si2) + hgfx(0)
    #     domain%m_fy(n0si2) = domain%m_fy(n0si2) + hgfy(0)
    #     domain%m_fz(n0si2) = domain%m_fz(n0si2) + hgfz(0)

    #     domain%m_fx(n1si2) = domain%m_fx(n1si2) + hgfx(1)
    #     domain%m_fy(n1si2) = domain%m_fy(n1si2) + hgfy(1)
    #     domain%m_fz(n1si2) = domain%m_fz(n1si2) + hgfz(1)

    #     domain%m_fx(n2si2) = domain%m_fx(n2si2) + hgfx(2)
    #     domain%m_fy(n2si2) = domain%m_fy(n2si2) + hgfy(2)
    #     domain%m_fz(n2si2) = domain%m_fz(n2si2) + hgfz(2)

    #     domain%m_fx(n3si2) = domain%m_fx(n3si2) + hgfx(3)
    #     domain%m_fy(n3si2) = domain%m_fy(n3si2) + hgfy(3)
    #     domain%m_fz(n3si2) = domain%m_fz(n3si2) + hgfz(3)

    #     domain%m_fx(n4si2) = domain%m_fx(n4si2) + hgfx(4)
    #     domain%m_fy(n4si2) = domain%m_fy(n4si2) + hgfy(4)
    #     domain%m_fz(n4si2) = domain%m_fz(n4si2) + hgfz(4)

    #     domain%m_fx(n5si2) = domain%m_fx(n5si2) + hgfx(5)
    #     domain%m_fy(n5si2) = domain%m_fy(n5si2) + hgfy(5)
    #     domain%m_fz(n5si2) = domain%m_fz(n5si2) + hgfz(5)

    #     domain%m_fx(n6si2) = domain%m_fx(n6si2) + hgfx(6)
    #     domain%m_fy(n6si2) = domain%m_fy(n6si2) + hgfy(6)
    #     domain%m_fz(n6si2) = domain%m_fz(n6si2) + hgfz(6)

    #     domain%m_fx(n7si2) = domain%m_fx(n7si2) + hgfx(7)
    #     domain%m_fy(n7si2) = domain%m_fy(n7si2) + hgfy(7)
    #     domain%m_fz(n7si2) = domain%m_fz(n7si2) + hgfz(7)
    # #endif
    end

    numNode = domain.numNode

    @inbounds for gnode in 1:numNode
        count = domain.nodeElemCount[gnode]
        start = domain.nodeElemStart[gnode]
        fx = 0.0
        fy = 0.0
        fz = 0.0
        @simd for i in 1:count
            elem = domain.nodeElemCornerList[start+i]
            fx += fx_elem[elem]
            fy += fy_elem[elem]
            fz += fz_elem[elem]
        end
        domain.fx[gnode] += fx
        domain.fy[gnode] += fy
        domain.fz[gnode] += fz
    end
end

function calcHourglassControlForElems(domain::Domain, determ, hgcoef)
    numElem = domain.numElem
    numElem8 = numElem * 8
    dvdx = Vector{Float64}(undef, numElem8)
    dvdy = Vector{Float64}(undef, numElem8)
    dvdz = Vector{Float64}(undef, numElem8)
    x8n = Vector{Float64}(undef, numElem8)
    y8n = Vector{Float64}(undef, numElem8)
    z8n = Vector{Float64}(undef, numElem8)

    #start loop over elements
    @inbounds for i in 1:numElem
        x1, y1, z1    = collectDomainNodesToElemNodes(domain, (i-1)*8+1)
        pfx, pfy, pfz = calcElemVolumeDerivative(x1, y1, z1)

        #   load into temporary storage for FB Hour Glass control
        for ii in 1:8
            jj = 8*(i-1) + ii

            dvdx[jj] = pfx[ii]
            dvdy[jj] = pfy[ii]
            dvdz[jj] = pfz[ii]
            x8n[jj]  = x1[ii]
            y8n[jj]  = y1[ii]
            z8n[jj]  = z1[ii]
        end

        determ[i] = domain.volo[i] * domain.v[i]

        #  Do a check for negative volumes
        if domain.v[i] <= 0.0
            error("Volume Error: Volume is negative")
        end
    end

    if hgcoef > 0.0
        MultiDcalcFBHourglassForceForElems(domain,determ,x8n,y8n,z8n,dvdx,dvdy,dvdz,hgcoef)
        # calcFBHourglassForceForElems(domain,determ,x8n,y8n,z8n,dvdx,dvdy,dvdz,hgcoef)
    end

    return nothing
end

function calcVolumeForceForElems(domain::Domain)
    # Based on FORTRAN implementation down from here
    hgcoef = domain.hgcoef
    numElem = domain.numElem
    VTD = typeof(domain.x)
    sigxx = VTD(undef, numElem)
    sigyy = VTD(undef, numElem)
    sigzz = VTD(undef, numElem)
    determ = VTD(undef, numElem)

    # Sum contributions to total stress tensor
    # initStressTermsForElems(domain, sigxx, sigyy, sigzz)

    #   call elemlib stress integration loop to produce nodal forces from
    #   material stresses.
    integrateStressForElems(domain, sigxx, sigyy, sigzz, determ)
    # MultiDintegrateStressForElems(domain, sigxx, sigyy, sigzz, determ)

    # check for negative element volume and abort if found
    for i in 1:numElem
        if determ[i] <= 0.0
            error("Mid Volume Error")
        end
    end
    calcHourglassControlForElems(domain, determ, hgcoef)
end

function calcForceForNodes(domain::Domain)
    commRecv(domain, MSG_COMM_SBN, 3,
             domain.sizeX + 1, domain.sizeY + 1, domain.sizeZ + 1,
             true, false)

    domain.fx .= 0.0
    domain.fy .= 0.0
    domain.fz .= 0.0

    calcVolumeForceForElems(domain);
    fields = (domain.fx, domain.fy, domain.fz)
    commSend(domain, MSG_COMM_SBN, fields,
             domain.sizeX + 1, domain.sizeY + 1, domain.sizeZ + 1,
             true, false)
    commSBN(domain, fields)
end

function calcAccelerationForNodes(domain::Domain)
    domain.xdd .= domain.fx ./ domain.nodalMass
    domain.ydd .= domain.fy ./ domain.nodalMass
    domain.zdd .= domain.fz ./ domain.nodalMass
end

function applyAccelerationBoundaryConditionsForNodes(domain::Domain)

    numNodeBC = (domain.sizeX+1)*(domain.sizeX+1)

    if length(domain.symmX) != 0
        for i in 1:numNodeBC
            domain.xdd[domain.symmX[i]+1] = 0.0
        end
    end
    if length(domain.symmY) != 0
        for i in 1:numNodeBC
            domain.ydd[domain.symmY[i]+1] = 0.0
        end
        end
    if length(domain.symmZ) != 0
        for i in 1:numNodeBC
            domain.zdd[domain.symmZ[i]+1] = 0.0
        end
    end
end

function calcVelocityForNodes(domain::Domain, dt, u_cut)


    numNode = domain.numNode

    for i in 1:numNode
        xdtmp = domain.xd[i] + domain.xdd[i] * dt
        if abs(xdtmp) < u_cut
            xdtmp = 0.0
        end
        domain.xd[i] = xdtmp

        ydtmp = domain.yd[i] + domain.ydd[i] * dt
        if abs(ydtmp) < u_cut
            ydtmp = 0.0
        end
        domain.yd[i] = ydtmp

        zdtmp = domain.zd[i] + domain.zdd[i] * dt
        if abs(zdtmp) < u_cut
            zdtmp = 0.0
        end
        domain.zd[i] = zdtmp
    end
end

function calcPositionForNodes(domain::Domain, dt)
    domain.x .= domain.x .+ domain.xd .* dt
    domain.y .= domain.y .+ domain.yd .* dt
    domain.z .= domain.z .+ domain.zd .* dt
end

function lagrangeNodal(domain::Domain)
    delt = domain.deltatime

    u_cut = domain.u_cut
    # time of boundary condition evaluation is beginning of step for force and
    # acceleration boundary conditions.
    calcForceForNodes(domain)

    if SEDOV_SYNC_POS_VEL_EARLY
        commRecv(domain, MSG_SYNC_POS_VEL, 6,
                 domain.sizeX + 1, domain.sizeY + 1, domain.sizeZ + 1,
                 false, false)
    end

    calcAccelerationForNodes(domain)

    applyAccelerationBoundaryConditionsForNodes(domain)

    calcVelocityForNodes(domain, delt, u_cut)
    calcPositionForNodes(domain, delt)

    if SEDOV_SYNC_POS_VEL_EARLY
        fields = (domain.x, domain.y, domain.z, domain.xd, domain.yd, domain.zd)
        commSend(domain, MSG_SYNC_POS_VEL, fields,
                 domain.sizeX + 1, domain.sizeY + 1, domain.sizeZ + 1,
                 false, false)
        # printAllFields(domain, "$(@__FILE__):$(@__LINE__)")
        commSyncPosVel(domain)
        # printAllFields(domain, "$(@__FILE__):$(@__LINE__)")
    end

    return nothing
end

@inline function areaFace( x0, x1, x2, x3,
                   y0, y1, y2, y3,
                   z0, z1, z2, z3  )

  fx = (x2 - x0) - (x3 - x1)
  fy = (y2 - y0) - (y3 - y1)
  fz = (z2 - z0) - (z3 - z1)
  gx = (x2 - x0) + (x3 - x1)
  gy = (y2 - y0) + (y3 - y1)
  gz = (z2 - z0) + (z3 - z1)

  area = (fx * fx + fy * fy + fz * fz) *
    (gx * gx + gy * gy + gz * gz) -
    (fx * gx + fy * gy + fz * gz) *
    (fx * gx + fy * gy + fz * gz)

  return area
end


@inline function JcalcElemCharacteristicLength( x, y, z, volume, i)

    charLength = 0.0

    a = areaFace(x[i + 1],x[i + 2],x[i + 3],x[i + 4],
                y[i + 1],y[i + 2],y[i + 3],y[i + 4],
                z[i + 1],z[i + 2],z[i + 3],z[i + 4])
    charLength = max(a,charLength)

    a = areaFace(x[i + 5],x[i + 6],x[i + 7],x[i + 8],
                y[i + 5],y[i + 6],y[i + 7],y[i + 8],
                z[i + 5],z[i + 6],z[i + 7],z[i + 8])
    charLength = max(a,charLength)

    a = areaFace(x[i + 1],x[i + 2],x[i + 6],x[i + 5],
                y[i + 1],y[i + 2],y[i + 6],y[i + 5],
                z[i + 1],z[i + 2],z[i + 6],z[i + 5])
    charLength = max(a,charLength)

    a = areaFace(x[i + 2],x[i + 3],x[i + 7],x[i + 6],
                y[i + 2],y[i + 3],y[i + 7],y[i + 6],
                z[i + 2],z[i + 3],z[i + 7],z[i + 6])
    charLength = max(a,charLength)

    a = areaFace(x[i + 3],x[i + 4],x[i + 8],x[i + 7],
                y[i + 3],y[i + 4],y[i + 8],y[i + 7],
                z[i + 3],z[i + 4],z[i + 8],z[i + 7])
    charLength = max(a,charLength)

    a = areaFace(x[i + 4],x[i + 1],x[i + 5],x[i + 8],
                 y[i + 4],y[i + 1],y[i + 5],y[i + 8],
                 z[i + 4],z[i + 1],z[i + 5],z[i + 8])
    charLength = max(a,charLength)

    charLength = 4.0 * volume / sqrt(charLength)

    return charLength
end

@inline function calcElemCharacteristicLength( x, y, z, volume)

    charLength = 0.0

    a = areaFace(x[1],x[2],x[3],x[4],
                y[1],y[2],y[3],y[4],
                z[1],z[2],z[3],z[4])
    charLength = max(a,charLength)

    a = areaFace(x[5],x[6],x[7],x[8],
                y[5],y[6],y[7],y[8],
                z[5],z[6],z[7],z[8])
    charLength = max(a,charLength)

    a = areaFace(x[1],x[2],x[6],x[5],
                y[1],y[2],y[6],y[5],
                z[1],z[2],z[6],z[5])
    charLength = max(a,charLength)

    a = areaFace(x[2],x[3],x[7],x[6],
                y[2],y[3],y[7],y[6],
                z[2],z[3],z[7],z[6])
    charLength = max(a,charLength)

    a = areaFace(x[3],x[4],x[8],x[7],
                y[3],y[4],y[8],y[7],
                z[3],z[4],z[8],z[7])
    charLength = max(a,charLength)

    a = areaFace(x[4],x[1],x[5],x[8],
                y[4],y[1],y[5],y[8],
                z[4],z[1],z[5],z[8])
    charLength = max(a,charLength)

    charLength = 4.0 * volume / sqrt(charLength)

    return charLength
end

@inline function calcElemShapeFunctionDerivatives(x, y, z)
    @inbounds begin
    x0 = x[1]
    x1 = x[2]
    x2 = x[3]
    x3 = x[4]
    x4 = x[5]
    x5 = x[6]
    x6 = x[7]
    x7 = x[8]

    y0 = y[1]
    y1 = y[2]
    y2 = y[3]
    y3 = y[4]
    y4 = y[5]
    y5 = y[6]
    y6 = y[7]
    y7 = y[8]

    z0 = z[1]
    z1 = z[2]
    z2 = z[3]
    z3 = z[4]
    z4 = z[5]
    z5 = z[6]
    z6 = z[7]
    z7 = z[8]

    fjxxi = .125 * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) )
    fjxet = .125 * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) )
    fjxze = .125 * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) )

    fjyxi = .125 * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) )
    fjyet = .125 * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) )
    fjyze = .125 * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) )

    fjzxi = .125 * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) )
    fjzet = .125 * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) )
    fjzze = .125 * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) )

    # compute cofactors
    cjxxi =    (fjyet * fjzze) - (fjzet * fjyze)
    cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze)
    cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet)

    cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze)
    cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze)
    cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet)

    cjzxi =    (fjxet * fjyze) - (fjyet * fjxze)
    cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze)
    cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet)

    # calculate partials :
    #     this need only be done for l = 0,1,2,3   since , by symmetry ,
    #     (6,7,4,5) = - (0,1,2,3) .
    b = MMatrix{8, 3, Float64}(undef) # shape function derivatives
    b[1,1] =   -  cjxxi  -  cjxet  -  cjxze
    b[2,1] =      cjxxi  -  cjxet  -  cjxze
    b[3,1] =      cjxxi  +  cjxet  -  cjxze
    b[4,1] =   -  cjxxi  +  cjxet  -  cjxze
    b[5,1] = -b[3,1]
    b[6,1] = -b[4,1]
    b[7,1] = -b[1,1]
    b[8,1] = -b[2,1]

    b[1,2] =   -  cjyxi  -  cjyet  -  cjyze
    b[2,2] =      cjyxi  -  cjyet  -  cjyze
    b[3,2] =      cjyxi  +  cjyet  -  cjyze
    b[4,2] =   -  cjyxi  +  cjyet  -  cjyze
    b[5,2] = -b[3,2]
    b[6,2] = -b[4,2]
    b[7,2] = -b[1,2]
    b[8,2] = -b[2,2]

    b[1,3] =   -  cjzxi  -  cjzet  -  cjzze
    b[2,3] =      cjzxi  -  cjzet  -  cjzze
    b[3,3] =      cjzxi  +  cjzet  -  cjzze
    b[4,3] =   -  cjzxi  +  cjzet  -  cjzze
    b[5,3] = -b[3,3]
    b[6,3] = -b[4,3]
    b[7,3] = -b[1,3]
    b[8,3] = -b[2,3]
    end #inbounds

    # calculate jacobian determinant (volume)
    el_volume = 8.0 * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet)
    return SMatrix(b), el_volume
end

@inline function calcElemVelocityGradient( xvel, yvel, zvel, b, detJ)
    @inbounds begin
    inv_detJ = 1.0 / detJ
    pfx = b[:,1]
    pfy = b[:,2]
    pfz = b[:,3]

    d1 = inv_detJ * (   pfx[1] * (xvel[1]-xvel[7])
                        + pfx[2] * (xvel[2]-xvel[8])
                        + pfx[3] * (xvel[3]-xvel[5])
                        + pfx[4] * (xvel[4]-xvel[6]) )

    d2 = inv_detJ * (   pfy[1] * (yvel[1]-yvel[7])
                        + pfy[2] * (yvel[2]-yvel[8])
                        + pfy[3] * (yvel[3]-yvel[5])
                        + pfy[4] * (yvel[4]-yvel[6]) )

    d3 = inv_detJ * (   pfz[1] * (zvel[1]-zvel[7])
                        + pfz[2] * (zvel[2]-zvel[8])
                        + pfz[3] * (zvel[3]-zvel[5])
                        + pfz[4] * (zvel[4]-zvel[6]) )

    dyddx = inv_detJ * (  pfx[1] * (yvel[1]-yvel[7])
                        + pfx[2] * (yvel[2]-yvel[8])
                        + pfx[3] * (yvel[3]-yvel[5])
                        + pfx[4] * (yvel[4]-yvel[6]) )

    dxddy = inv_detJ * (  pfy[1] * (xvel[1]-xvel[7])
                        + pfy[2] * (xvel[2]-xvel[8])
                        + pfy[3] * (xvel[3]-xvel[5])
                        + pfy[4] * (xvel[4]-xvel[6]) )

    dzddx = inv_detJ * (  pfx[1] * (zvel[1]-zvel[7])
                        + pfx[2] * (zvel[2]-zvel[8])
                        + pfx[3] * (zvel[3]-zvel[5])
                        + pfx[4] * (zvel[4]-zvel[6]) )

    dxddz = inv_detJ * (  pfz[1] * (xvel[1]-xvel[7])
                        + pfz[2] * (xvel[2]-xvel[8])
                        + pfz[3] * (xvel[3]-xvel[5])
                        + pfz[4] * (xvel[4]-xvel[6]) )

    dzddy = inv_detJ * (  pfy[1] * (zvel[1]-zvel[7])
                        + pfy[2] * (zvel[2]-zvel[8])
                        + pfy[3] * (zvel[3]-zvel[5])
                        + pfy[4] * (zvel[4]-zvel[6]) )

    dyddz = inv_detJ * (  pfz[1] * (yvel[1]-yvel[7])
                        + pfz[2] * (yvel[2]-yvel[8])
                        + pfz[3] * (yvel[3]-yvel[5])
                        + pfz[4] * (yvel[4]-yvel[6]) )
    end #inbounds

    d6 = 0.5 * ( dxddy + dyddx )
    d5 = 0.5 * ( dxddz + dzddx )
    d4 = 0.5 * ( dzddy + dyddz )

    return SVector(d1, d2, d3, d4, d5, d6)
end

@inline function collectNodal(nodelist, src, i)
    @inbounds begin
        s1 = src[nodelist[i+1]+1]
        s2 = src[nodelist[i+2]+1]
        s3 = src[nodelist[i+3]+1]
        s4 = src[nodelist[i+4]+1]
        s5 = src[nodelist[i+5]+1]
        s6 = src[nodelist[i+6]+1]
        s7 = src[nodelist[i+7]+1]
        s8 = src[nodelist[i+8]+1]
    end

    return SVector(s1, s2, s3, s4, s5, s6, s7, s8)
end

function calcKinematicsForElems_A(i, x_single, y_single, z_single, volo, v, vnew, delv, arealg)

    i8 = 8 * (i - 1)
    volume = JcalcElemVolume(x_single, y_single, z_single, i8)
    relativeVolume = volume / volo[i]
    vnew[i] = relativeVolume

    delv[i] = relativeVolume - v[i]
    
    arealg[i] = JcalcElemCharacteristicLength(x_single, y_single, z_single, volume, i8)

end

function calcKinematicsForElems_B()

    i8 = 8 * (i - 1)

    dt2 = 0.5 * dt

    for j in 1:8
        x_single[i8 + j] = x_single[i8 + j] - dt2 * xd_single[i8 + j]
        y_single[i8 + j] = y_single[i8 + j] - dt2 * yd_single[i8 + j]
        z_single[i8 + j] = z_single[i8 + j] - dt2 * zd_single[i8 + j]
    end

    # for i in 1:8:length(x_single)
    #     x_single[i] -= dt2 * xd_single[i]
    #     y_single[i] -= dt2 * yd_single[i]
    #     z_single[i] -= dt2 * zd_single[i]
    # end
end

function calcKinematicsForElems_C(i, dxx, dyy, dzz, D)

    i6 = 6 * (i - 1)

    dxx[i] = D[i6 + 1]
    dyy[i] = D[i6 + 2]
    dzz[i] = D[i6 + 3]

end

@inline function JcalcElemVelocityGradient(i, xvel, yvel, zvel, b, detJ, D)
    @inbounds begin

        i24 = 24 * (i - 1)
        i8 = 8 * (i - 1) 
        i6 = 6 * (i - 1)

        if detJ[i] <= 0.0
            error("Element has zero or negative volume")
        end

        inv_detJ = 1.0 / detJ[i]
        

        d1 = inv_detJ * (       b[i24 + 1] * ( xvel[i8 + 1] - xvel[i8+7])
                            +   b[i24 + 2] * ( xvel[i8 + 2] - xvel[i8+8])
                            +   b[i24 + 3] * ( xvel[i8 + 3] - xvel[i8+5])
                            +   b[i24 + 4] * ( xvel[i8 + 4] - xvel[i8+6]) )

        d2 = inv_detJ * (       b[i24 + 9]  * ( yvel[i8 + 1] - yvel[i8+7])
                            +   b[i24 + 10] * ( yvel[i8 + 2] - yvel[i8+8])
                            +   b[i24 + 11] * ( yvel[i8 + 3] - yvel[i8+5])
                            +   b[i24 + 12] * ( yvel[i8 + 4] - yvel[i8+6]) )

        d3 = inv_detJ * (       b[i24 + 17] * ( zvel[i8 + 1] - zvel[i8+7])
                            +   b[i24 + 18] * ( zvel[i8 + 2] - zvel[i8+8])
                            +   b[i24 + 19] * ( zvel[i8 + 3] - zvel[i8+5])
                            +   b[i24 + 20] * ( zvel[i8 + 4] - zvel[i8+6]) )

        dyddx = inv_detJ * (    b[i24 + 1] * ( yvel[i8 + 1] - yvel[i8+7])
                            +   b[i24 + 2] * ( yvel[i8 + 2] - yvel[i8+8])
                            +   b[i24 + 3] * ( yvel[i8 + 3] - yvel[i8+5])
                            +   b[i24 + 4] * ( yvel[i8 + 4] - yvel[i8+6]) )

        dxddy = inv_detJ * (    b[i24 + 9]  * ( xvel[i8 + 1] - xvel[i8+7])
                            +   b[i24 + 10] * ( xvel[i8 + 2] - xvel[i8+8])
                            +   b[i24 + 11] * ( xvel[i8 + 3] - xvel[i8+5])
                            +   b[i24 + 12] * ( xvel[i8 + 4] - xvel[i8+6]) )

        dzddx = inv_detJ * (    b[i24 + 1] * ( zvel[i8 + 1] - zvel[i8+7])
                            +   b[i24 + 2] * ( zvel[i8 + 2] - zvel[i8+8])
                            +   b[i24 + 3] * ( zvel[i8 + 3] - zvel[i8+5])
                            +   b[i24 + 4] * ( zvel[i8 + 4] - zvel[i8+6]) )

        dxddz = inv_detJ * (    b[i24 + 17] * ( xvel[i8 + 1] - xvel[i8+7])
                            +   b[i24 + 18] * ( xvel[i8 + 2] - xvel[i8+8])
                            +   b[i24 + 19] * ( xvel[i8 + 3] - xvel[i8+5])
                            +   b[i24 + 20] * ( xvel[i8 + 4] - xvel[i8+6]) )

        dzddy = inv_detJ * (    b[i24 + 9] *  ( zvel[i8 + 1] - zvel[i8+7])
                            +   b[i24 + 10] * ( zvel[i8 + 2] - zvel[i8+8])
                            +   b[i24 + 11] * ( zvel[i8 + 3] - zvel[i8+5])
                            +   b[i24 + 12] * ( zvel[i8 + 4] - zvel[i8+6]) )

        dyddz = inv_detJ * (    b[i24 + 17] * ( yvel[i8 + 1] - yvel[i8+7])
                            +   b[i24 + 18] * ( yvel[i8 + 2] - yvel[i8+8])
                            +   b[i24 + 19] * ( yvel[i8 + 3] - yvel[i8+5])
                            +   b[i24 + 20] * ( yvel[i8 + 4] - yvel[i8+6]) )

    end #inbounds
    
    
    d6 = 0.5 * ( dxddy + dyddx )
    d5 = 0.5 * ( dxddz + dzddx )
    d4 = 0.5 * ( dzddy + dyddz )

    D[i6 + 1] = d1
    D[i6 + 2] = d2
    D[i6 + 3] = d3
    D[i6 + 4] = d4
    D[i6 + 5] = d5
    D[i6 + 6] = d6

    #  return SVector(d1, d2, d3, d4, d5, d6)
end

@inline function calcKinematicsForElems_Aa(i, dt, x_single, y_single, z_single, xd_single, yd_single, zd_single)

    i8 = 8 * (i - 1)
    dt2 = 0.5 * dt
    for j in 1:8
        x_single[i8 + j] -= dt2 * xd_single[i8 + j]
        y_single[i8 + j] -= dt2 * yd_single[i8 + j]
        z_single[i8 + j] -= dt2 * zd_single[i8 + j]
    end

end

function MultiDcalcKinematicsForElems(domain::Domain, numElem, dt)


    numElem8 = numElem * 8
    numElem3 = numElem * 3

    # Initial Arrays
    x_single = Vector{Float64}(undef, numElem8)
    y_single = Vector{Float64}(undef, numElem8)
    z_single = Vector{Float64}(undef, numElem8)
    xd_single = Vector{Float64}(undef, numElem8)
    yd_single = Vector{Float64}(undef, numElem8)
    zd_single = Vector{Float64}(undef, numElem8)
    b_single = Vector{Float64}(undef, numElem3 * 8)
    determ = Vector{Float64}(undef, numElem)
    D = Vector{Float64}(undef, numElem * 6)
    # Dxx = Vector{Float64}(undef, numElem)
    # Dyy = Vector{Float64}(undef, numElem)
    # Dzz = Vector{Float64}(undef, numElem)

    # JACC Arrays
    nodelist_d = JACC.Array(domain.nodelist)
    x_d = JACC.Array(domain.x)
    y_d = JACC.Array(domain.y)
    z_d = JACC.Array(domain.z)
    xd_d = JACC.Array(domain.xd)
    yd_d = JACC.Array(domain.yd)
    zd_d = JACC.Array(domain.zd)
    volo_d = JACC.Array(domain.volo)
    vnew_d = JACC.Array(domain.vnew)
    determ_d = JACC.Array(determ)
    v_d = JACC.Array(domain.v)
    delv_d = JACC.Array(domain.delv)
    arealg_d = JACC.Array(domain.arealg)
    dxx_d = JACC.Array(domain.dxx)
    dyy_d = JACC.Array(domain.dyy)
    dzz_d = JACC.Array(domain.dzz)
    # Dxx_d = JACC.Array(Dxx)
    # Dyy_d = JACC.Array(Dyy)
    # Dzz_d = JACC.Array(Dzz)
    

    # JACC Big Arrays
    x_single_d = JACC.Array(x_single)
    y_single_d = JACC.Array(y_single)
    z_single_d = JACC.Array(z_single)
    xd_single_d = JACC.Array(xd_single)
    yd_single_d = JACC.Array(yd_single)
    zd_single_d = JACC.Array(zd_single)
    b_single_d = JACC.Array(b_single)
    D_d = JACC.Array(D)

    time1 = time()
    JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, x_single_d, x_d)
    JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, y_single_d, y_d)
    JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, z_single_d, z_d)
    JACC.parallel_for(domain.numElem, calcKinematicsForElems_A, x_single_d, y_single_d, z_single_d, volo_d, v_d, vnew_d, delv_d, arealg_d)
    JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, xd_single_d, xd_d)
    JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, yd_single_d, yd_d)
    JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, zd_single_d, zd_d)
    JACC.parallel_for(domain.numElem, calcKinematicsForElems_Aa, dt, x_single_d, y_single_d, z_single_d, xd_single_d, yd_single_d, zd_single_d)
    JACC.parallel_for(domain.numElem, JcalcElemShapeFunctionDerivatives, x_single_d, y_single_d, z_single_d, b_single_d, determ_d) 
    JACC.parallel_for(domain.numElem, JcalcElemVelocityGradient, xd_single_d, yd_single_d, zd_single_d, b_single_d, determ_d, D_d)
    JACC.parallel_for(domain.numElem, calcKinematicsForElems_C, dxx_d, dyy_d, dzz_d, D_d)
    time2 = time()
    println("Time Taken for MultiDCalcKinematicsForElems: ", time2 - time1)
    
    
    total_time_collectNodal_x = 0.0
    total_time_collectNodal_y = 0.0
    total_time_collectNodal_z = 0.0
    total_time_calcKinematicsForElems_A = 0.0
    total_time_collectNodal_xd = 0.0
    total_time_collectNodal_yd = 0.0
    total_time_collectNodal_zd = 0.0
    total_time_calcKinematicsForElems_Aa = 0.0
    total_time_calcElemShapeFunctionDerivatives = 0.0
    total_time_calcElemVelocityGradient = 0.0
    total_time_calcKinematicsForElems_C = 0.0

    num_runs = 10

    for i in 1:num_runs
        time1 = time()
        JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, x_single_d, x_d)
        time2 = time()
        total_time_collectNodal_x += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, y_single_d, y_d)
        time2 = time()
        total_time_collectNodal_y += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, z_single_d, z_d)
        time2 = time()
        total_time_collectNodal_z += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, calcKinematicsForElems_A, x_single_d, y_single_d, z_single_d, volo_d, v_d, vnew_d, delv_d, arealg_d)
        time2 = time()
        total_time_calcKinematicsForElems_A += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, xd_single_d, xd_d)
        time2 = time()
        total_time_collectNodal_xd += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, yd_single_d, yd_d)
        time2 = time()
        total_time_collectNodal_yd += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, JcollectNodal, nodelist_d, zd_single_d, zd_d)
        time2 = time()
        total_time_collectNodal_zd += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, calcKinematicsForElems_Aa, dt, x_single_d, y_single_d, z_single_d, xd_single_d, yd_single_d, zd_single_d)
        time2 = time()
        total_time_calcKinematicsForElems_Aa += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, JcalcElemShapeFunctionDerivatives, x_single_d, y_single_d, z_single_d, b_single_d, determ_d)
        time2 = time()
        total_time_calcElemShapeFunctionDerivatives += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, JcalcElemVelocityGradient, xd_single_d, yd_single_d, zd_single_d, b_single_d, determ_d, D_d)
        time2 = time()
        total_time_calcElemVelocityGradient += (time2 - time1)

        time1 = time()
        JACC.parallel_for(domain.numElem, calcKinematicsForElems_C, dxx_d, dyy_d, dzz_d, D_d)
        time2 = time()
        total_time_calcKinematicsForElems_C += (time2 - time1)
    end

    average_time_collectNodal_x = total_time_collectNodal_x / num_runs
    average_time_collectNodal_y = total_time_collectNodal_y / num_runs
    average_time_collectNodal_z = total_time_collectNodal_z / num_runs
    average_time_calcKinematicsForElems_A = total_time_calcKinematicsForElems_A / num_runs
    average_time_collectNodal_xd = total_time_collectNodal_xd / num_runs
    average_time_collectNodal_yd = total_time_collectNodal_yd / num_runs
    average_time_collectNodal_zd = total_time_collectNodal_zd / num_runs
    average_time_calcKinematicsForElems_Aa = total_time_calcKinematicsForElems_Aa / num_runs
    average_time_calcElemShapeFunctionDerivatives = total_time_calcElemShapeFunctionDerivatives / num_runs
    average_time_calcElemVelocityGradient = total_time_calcElemVelocityGradient / num_runs
    average_time_calcKinematicsForElems_C = total_time_calcKinematicsForElems_C / num_runs

    total_average_time = average_time_collectNodal_x + average_time_collectNodal_y + average_time_collectNodal_z + average_time_calcKinematicsForElems_A + average_time_collectNodal_xd + average_time_collectNodal_yd + average_time_collectNodal_zd + average_time_calcKinematicsForElems_Aa + average_time_calcElemShapeFunctionDerivatives + average_time_calcElemVelocityGradient + average_time_calcKinematicsForElems_C

    println("Average time taken for JcollectNodal (x) over $num_runs runs: ", average_time_collectNodal_x)
    println("Average time taken for JcollectNodal (y) over $num_runs runs: ", average_time_collectNodal_y)
    println("Average time taken for JcollectNodal (z) over $num_runs runs: ", average_time_collectNodal_z)
    println("Average time taken for calcKinematicsForElems_A over $num_runs runs: ", average_time_calcKinematicsForElems_A)
    println("Average time taken for JcollectNodal (xd) over $num_runs runs: ", average_time_collectNodal_xd)
    println("Average time taken for JcollectNodal (yd) over $num_runs runs: ", average_time_collectNodal_yd)
    println("Average time taken for JcollectNodal (zd) over $num_runs runs: ", average_time_collectNodal_zd)
    println("Average time taken for calcKinematicsForElems_Aa over $num_runs runs: ", average_time_calcKinematicsForElems_Aa)
    println("Average time taken for JcalcElemShapeFunctionDerivatives over $num_runs runs: ", average_time_calcElemShapeFunctionDerivatives)
    println("Average time taken for JcalcElemVelocityGradient over $num_runs runs: ", average_time_calcElemVelocityGradient)
    println("Average time taken for calcKinematicsForElems_C over $num_runs runs: ", average_time_calcKinematicsForElems_C)
    println("Total average time taken for all functions over $num_runs runs: ", total_average_time)
    
    copyto!(domain.dxx, dxx_d)
    copyto!(domain.dyy, dyy_d)
    copyto!(domain.dzz, dzz_d)
    copyto!(domain.vnew, vnew_d)
    copyto!(domain.delv, delv_d)
    copyto!(domain.arealg, arealg_d)

end


function calcKinematicsForElems(domain::Domain, numElem, dt)

    nodelist = domain.nodelist
    # loop over all elements
    # printAllFields(domain, "$(@__FILE__):$(@__LINE__)")
    for k in 1:numElem
        # get nodal coordinates from global arrays and copy into local arrays
        @inbounds begin
            x_local = collectNodal(nodelist, domain.x, (k-1)*8)
            y_local = collectNodal(nodelist, domain.y, (k-1)*8)
            z_local = collectNodal(nodelist, domain.z, (k-1)*8)
        end

        # volume calculations
        volume = calcElemVolume(x_local, y_local, z_local)
        relativeVolume = volume / domain.volo[k]
        domain.vnew[k] = relativeVolume
        if domain.vnew[k] <= 0.0
            # @error "negative volume found in calcKinematicsForElems" rank=getMyRank(domain.comm) k volume volo=domain.volo[k] xyz=(x_local, y_local, z_local)
            error("Volume Error: calcKinematicsForElems")
        end
        domain.delv[k] = relativeVolume - domain.v[k]

        # set characteristic length
        domain.arealg[k] = calcElemCharacteristicLength(x_local, y_local,
                                                        z_local, volume)

        #  get nodal velocities from global array and copy into local arrays.
        @inbounds begin
            xd_local = collectNodal(nodelist, domain.xd, (k-1)*8)
            yd_local = collectNodal(nodelist, domain.yd, (k-1)*8)
            zd_local = collectNodal(nodelist, domain.zd, (k-1)*8)
        end

        dt2 = 0.5 * dt
        x_local = x_local .- dt2 .* xd_local
        y_local = y_local .- dt2 .* yd_local
        z_local = z_local .- dt2 .* zd_local

        B, detJ = calcElemShapeFunctionDerivatives(x_local, y_local, z_local)

        D = calcElemVelocityGradient(xd_local, yd_local, zd_local, B, detJ)

        # put velocity gradient quantities into their global arrays.
        domain.dxx[k] = D[1]
        domain.dyy[k] = D[2]
        domain.dzz[k] = D[3]
    end
end

function calcLagrangeElements(domain, delt)
    numElem = domain.numElem
    if numElem > 0
        # calcKinematicsForElems(domain, numElem, delt)
        MultiDcalcKinematicsForElems(domain, numElem, delt)

        # element loop to do some stuff not included in the elemlib function.

        for k in 1:numElem
        # calc strain rate and apply as constraint (only done in FB element)
            vdov = domain.dxx[k] + domain.dyy[k] + domain.dzz[k]
            vdovthird = vdov/3.0

            # make the rate of deformation tensor deviatoric
            domain.vdov[k] = vdov
            domain.dxx[k] = domain.dxx[k] - vdovthird
            domain.dyy[k] = domain.dyy[k] - vdovthird
            domain.dzz[k] = domain.dzz[k] - vdovthird

            # See if any volumes are negative, and take appropriate action.
            if domain.vnew[k] <= 0.0
                error("Volume Error :2157")
            end
        end
    end
end

function MultiDcalcMonotonicQGradientsForElems(i, nodelist, x_single, y_single, z_single, x, y, z, xv, yv, zv, xd, yd, zd, volo, vnew, delx_zeta, delv_zeta, delx_xi, delv_xi, delx_eta, delv_eta, ptiny, zero_array)  

    # zero_shared = JACC.experimental.shared(zero_array)
    # x_shared = JACC.experimental.shared(x)
    k = (i - 1) * 8



    n0 = nodelist[k + 1]
    n1 = nodelist[k + 2]
    n2 = nodelist[k + 3]
    n3 = nodelist[k + 4]
    n4 = nodelist[k + 5]
    n5 = nodelist[k + 6]
    n6 = nodelist[k + 7]
    n7 = nodelist[k + 8]

    x_single[k + 1] = x[n0 + 1]
    x_single[k + 2] = x[n1 + 1]
    x_single[k + 3] = x[n2 + 1]
    x_single[k + 4] = x[n3 + 1]
    x_single[k + 5] = x[n4 + 1]
    x_single[k + 6] = x[n5 + 1]
    x_single[k + 7] = x[n6 + 1]
    x_single[k + 8] = x[n7 + 1]

    # x_single0 = x[n0 + 1]
    # x_single1 = x[n1 + 1]
    # x_single2 = x[n2 + 1]
    # x_single3 = x[n3 + 1]
    # x_single4 = x[n4 + 1]
    # x_single5 = x[n5 + 1]
    # x_single6 = x[n6 + 1]
    # x_single7 = x[n7 + 1]

    y_single[k + 1] = y[n0 + 1]
    y_single[k + 2] = y[n1 + 1]
    y_single[k + 3] = y[n2 + 1]
    y_single[k + 4] = y[n3 + 1]
    y_single[k + 5] = y[n4 + 1]
    y_single[k + 6] = y[n5 + 1]
    y_single[k + 7] = y[n6 + 1]
    y_single[k + 8] = y[n7 + 1]

    # y_single0 = y[n0 + 1]
    # y_single1 = y[n1 + 1]
    # y_single2 = y[n2 + 1]
    # y_single3 = y[n3 + 1]
    # y_single4 = y[n4 + 1]
    # y_single5 = y[n5 + 1]
    # y_single6 = y[n6 + 1]
    # y_single7 = y[n7 + 1]

    z_single[k + 1] = z[n0 + 1]
    z_single[k + 2] = z[n1 + 1]
    z_single[k + 3] = z[n2 + 1]
    z_single[k + 4] = z[n3 + 1]
    z_single[k + 5] = z[n4 + 1]
    z_single[k + 6] = z[n5 + 1]
    z_single[k + 7] = z[n6 + 1]
    z_single[k + 8] = z[n7 + 1]

    # z_single0 = z[n0 + 1]
    # z_single1 = z[n1 + 1]
    # z_single2 = z[n2 + 1]
    # z_single3 = z[n3 + 1]
    # z_single4 = z[n4 + 1]
    # z_single5 = z[n5 + 1]
    # z_single6 = z[n6 + 1]
    # z_single7 = z[n7 + 1]

    xv[k + 1] = xd[n0 + 1]
    xv[k + 2] = xd[n1 + 1]
    xv[k + 3] = xd[n2 + 1]
    xv[k + 4] = xd[n3 + 1]
    xv[k + 5] = xd[n4 + 1]
    xv[k + 6] = xd[n5 + 1]
    xv[k + 7] = xd[n6 + 1]
    xv[k + 8] = xd[n7 + 1]

    # xv0 = xd[n0 + 1]
    # xv1 = xd[n1 + 1]
    # xv2 = xd[n2 + 1]
    # xv3 = xd[n3 + 1]
    # xv4 = xd[n4 + 1]
    # xv5 = xd[n5 + 1]
    # xv6 = xd[n6 + 1]
    # xv7 = xd[n7 + 1]

    yv[k + 1] = yd[n0 + 1]
    yv[k + 2] = yd[n1 + 1]
    yv[k + 3] = yd[n2 + 1]
    yv[k + 4] = yd[n3 + 1]
    yv[k + 5] = yd[n4 + 1]
    yv[k + 6] = yd[n5 + 1]
    yv[k + 7] = yd[n6 + 1]
    yv[k + 8] = yd[n7 + 1]

    # yv0 = yd[n0 + 1]
    # yv1 = yd[n1 + 1]
    # yv2 = yd[n2 + 1]
    # yv3 = yd[n3 + 1]
    # yv4 = yd[n4 + 1]
    # yv5 = yd[n5 + 1]
    # yv6 = yd[n6 + 1]
    # yv7 = yd[n7 + 1]

    zv[k + 1] = zd[n0 + 1]
    zv[k + 2] = zd[n1 + 1]
    zv[k + 3] = zd[n2 + 1]
    zv[k + 4] = zd[n3 + 1]
    zv[k + 5] = zd[n4 + 1]
    zv[k + 6] = zd[n5 + 1]
    zv[k + 7] = zd[n6 + 1]
    zv[k + 8] = zd[n7 + 1]

    # zv0 = zd[n0 + 1]
    # zv1 = zd[n1 + 1]
    # zv2 = zd[n2 + 1]
    # zv3 = zd[n3 + 1]
    # zv4 = zd[n4 + 1]
    # zv5 = zd[n5 + 1]
    # zv6 = zd[n6 + 1]
    # zv7 = zd[n7 + 1]

    vol = volo[i] * vnew[i]
    norm = 1.0 / (vol + ptiny)

    # Part 02 Begin



    dxj = -0.25 * (sum4(x_single[k + 1], x_single[k + 2], x_single[k + 6], x_single[k + 5]) - sum4(x_single[k + 4], x_single[k + 3], x_single[k + 7], x_single[k + 8]))
    dyj = -0.25 * (sum4(y_single[k + 1], y_single[k + 2], y_single[k + 6], y_single[k + 5]) - sum4(y_single[k + 4], y_single[k + 3], y_single[k + 7], y_single[k + 8]))
    dzj = -0.25 * (sum4(z_single[k + 1], z_single[k + 2], z_single[k + 6], z_single[k + 5]) - sum4(z_single[k + 4], z_single[k + 3], z_single[k + 7], z_single[k + 8]))

    # dxj = -0.25 * (sum4(x_single0, x_single1, x_single5, x_single4) - sum4(x_single3, x_single2, x_single6, x_single7))
    # dyj = -0.25 * (sum4(y_single0, y_single1, y_single5, y_single4) - sum4(y_single3, y_single2, y_single6, y_single7))
    # dzj = -0.25 * (sum4(z_single0, z_single1, z_single5, z_single4) - sum4(z_single3, z_single2, z_single6, z_single7))

    dxi = 0.25 * (sum4(x_single[k + 2], x_single[k + 3], x_single[k + 7], x_single[k + 6]) - sum4(x_single[k + 1], x_single[k + 4], x_single[k + 8], x_single[k + 5]))
    dyi = 0.25 * (sum4(y_single[k + 2], y_single[k + 3], y_single[k + 7], y_single[k + 6]) - sum4(y_single[k + 1], y_single[k + 4], y_single[k + 8], y_single[k + 5]))
    dzi = 0.25 * (sum4(z_single[k + 2], z_single[k + 3], z_single[k + 7], z_single[k + 6]) - sum4(z_single[k + 1], z_single[k + 4], z_single[k + 8], z_single[k + 5]))

    # dxi = 0.25 * (sum4(x_single1, x_single2, x_single6, x_single5) - sum4(x_single0, x_single3, x_single7, x_single4))
    # dyi = 0.25 * (sum4(y_single1, y_single2, y_single6, y_single5) - sum4(y_single0, y_single3, y_single7, y_single4))
    # dzi = 0.25 * (sum4(z_single1, z_single2, z_single6, z_single5) - sum4(z_single0, z_single3, z_single7, z_single4))

    dxk = 0.25 * (sum4(x_single[k + 5], x_single[k + 6], x_single[k + 7], x_single[k + 8]) - sum4(x_single[k + 1], x_single[k + 2], x_single[k + 3], x_single[k + 4]))
    dyk = 0.25 * (sum4(y_single[k + 5], y_single[k + 6], y_single[k + 7], y_single[k + 8]) - sum4(y_single[k + 1], y_single[k + 2], y_single[k + 3], y_single[k + 4]))
    dzk = 0.25 * (sum4(z_single[k + 5], z_single[k + 6], z_single[k + 7], z_single[k + 8]) - sum4(z_single[k + 1], z_single[k + 2], z_single[k + 3], z_single[k + 4]))

    # dxk = 0.25 * (sum4(x_single4, x_single5, x_single6, x_single7) - sum4(x_single0, x_single1, x_single2, x_single3))
    # dyk = 0.25 * (sum4(y_single4, y_single5, y_single6, y_single7) - sum4(y_single0, y_single1, y_single2, y_single3))
    # dzk = 0.25 * (sum4(z_single4, z_single5, z_single6, z_single7) - sum4(z_single0, z_single1, z_single2, z_single3))



    ax = dyi * dzj - dzi * dyj
    ay = dzi * dxj - dxi * dzj
    az = dxi * dyj - dyi * dxj
    # HERE 01
    delx_zeta[i] = vol / sqrt(ax * ax + ay * ay + az * az + ptiny)

    ax = ax * norm
    ay = ay * norm
    az = az * norm

    dxv = 0.25 * (sum4(xv[k + 5], xv[k + 6], xv[k + 7], xv[k + 8]) - sum4(xv[k + 1], xv[k + 2], xv[k + 3], xv[k + 4]))
    dyv = 0.25 * (sum4(yv[k + 5], yv[k + 6], yv[k + 7], yv[k + 8]) - sum4(yv[k + 1], yv[k + 2], yv[k + 3], yv[k + 4]))
    dzv = 0.25 * (sum4(zv[k + 5], zv[k + 6], zv[k + 7], zv[k + 8]) - sum4(zv[k + 1], zv[k + 2], zv[k + 3], zv[k + 4]))

    # dxv = 0.25 * (sum4(xv4, xv5, xv6, xv7) - sum4(xv0, xv1, xv2, xv3))
    # dyv = 0.25 * (sum4(yv4, yv5, yv6, yv7) - sum4(yv0, yv1, yv2, yv3))
    # dzv = 0.25 * (sum4(zv4, zv5, zv6, zv7) - sum4(zv0, zv1, zv2, zv3))

    delv_zeta[i] = ax * dxv + ay * dyv + az * dzv 

    # Part 02 End
    
    # Part 03 Begin

    ax = dyj * dzk - dzj * dyk
    ay = dzj * dxk - dxj * dzk
    az = dxj * dyk - dyj * dxk

    delx_xi[i] = vol / sqrt(ax * ax + ay * ay + az * az + ptiny)

    ax = ax * norm
    ay = ay * norm
    az = az * norm
    # HERE 02

    dxv = 0.25 * (sum4(xv[k + 2], xv[k + 3], xv[k + 7], xv[k + 6]) - sum4(xv[k + 1], xv[k + 4], xv[k + 8], xv[k + 5]))
    dyv = 0.25 * (sum4(yv[k + 2], yv[k + 3], yv[k + 7], yv[k + 6]) - sum4(yv[k + 1], yv[k + 4], yv[k + 8], yv[k + 5]))
    dzv = 0.25 * (sum4(zv[k + 2], zv[k + 3], zv[k + 7], zv[k + 6]) - sum4(zv[k + 1], zv[k + 4], zv[k + 8], zv[k + 5]))

    # dxv = 0.25 * (sum4(xv1, xv2, xv6, xv5) - sum4(xv0, xv3, xv7, xv4))
    # dyv = 0.25 * (sum4(yv1, yv2, yv6, yv5) - sum4(yv0, yv3, yv7, yv4))
    # dzv = 0.25 * (sum4(zv1, zv2, zv6, zv5) - sum4(zv0, zv3, zv7, zv4))

    delv_xi[i] = ax * dxv + ay * dyv + az * dzv

    ax = dyk * dzi - dzk * dyi ;
    ay = dzk * dxi - dxk * dzi ; 
    az = dxk * dyi - dyk * dxi ;

    delx_eta[i] = vol / sqrt(ax * ax + ay * ay + az * az + ptiny)

    ax = ax * norm
    ay = ay * norm
    az = az * norm

    dxv = -0.25 * (sum4(xv[k + 1], xv[k + 2], xv[k + 6], xv[k + 5]) - sum4(xv[k + 4], xv[k + 3], xv[k + 7], xv[k + 8]))
    dyv = -0.25 * (sum4(yv[k + 1], yv[k + 2], yv[k + 6], yv[k + 5]) - sum4(yv[k + 4], yv[k + 3], yv[k + 7], yv[k + 8]))
    dzv = -0.25 * (sum4(zv[k + 1], zv[k + 2], zv[k + 6], zv[k + 5]) - sum4(zv[k + 4], zv[k + 3], zv[k + 7], zv[k + 8]))

    # dxv = -0.25 * (sum4(xv0, xv1, xv5, xv4) - sum4(xv3, xv2, xv6, xv7))
    # dyv = -0.25 * (sum4(yv0, yv1, yv5, yv4) - sum4(yv3, yv2, yv6, yv7))
    # dzv = -0.25 * (sum4(zv0, zv1, zv5, zv4) - sum4(zv3, zv2, zv6, zv7))

    delv_eta[i] = ax * dxv + ay * dyv + az * dzv

    # Part 03 End


end

@inline function sum4(x1, x2, x3, x4)
    x1 + x2 + x3 + x4
end

function calcMonotonicQGradientsForElems(domain::Domain)

  ptiny = 1.e-36

  numElem = domain.numElem

  for i in 1:numElem
        k = (i-1)*8
        n0 = domain.nodelist[k+1]
        n1 = domain.nodelist[k+2]
        n2 = domain.nodelist[k+3]
        n3 = domain.nodelist[k+4]
        n4 = domain.nodelist[k+5]
        n5 = domain.nodelist[k+6]
        n6 = domain.nodelist[k+7]
        n7 = domain.nodelist[k+8]

        x0 = domain.x[n0+1]
        x1 = domain.x[n1+1]
        x2 = domain.x[n2+1]
        x3 = domain.x[n3+1]
        x4 = domain.x[n4+1]
        x5 = domain.x[n5+1]
        x6 = domain.x[n6+1]
        x7 = domain.x[n7+1]

        y0 = domain.y[n0+1]
        y1 = domain.y[n1+1]
        y2 = domain.y[n2+1]
        y3 = domain.y[n3+1]
        y4 = domain.y[n4+1]
        y5 = domain.y[n5+1]
        y6 = domain.y[n6+1]
        y7 = domain.y[n7+1]

        z0 = domain.z[n0+1]
        z1 = domain.z[n1+1]
        z2 = domain.z[n2+1]
        z3 = domain.z[n3+1]
        z4 = domain.z[n4+1]
        z5 = domain.z[n5+1]
        z6 = domain.z[n6+1]
        z7 = domain.z[n7+1]

        xv0 = domain.xd[n0+1]
        xv1 = domain.xd[n1+1]
        xv2 = domain.xd[n2+1]
        xv3 = domain.xd[n3+1]
        xv4 = domain.xd[n4+1]
        xv5 = domain.xd[n5+1]
        xv6 = domain.xd[n6+1]
        xv7 = domain.xd[n7+1]

        yv0 = domain.yd[n0+1]
        yv1 = domain.yd[n1+1]
        yv2 = domain.yd[n2+1]
        yv3 = domain.yd[n3+1]
        yv4 = domain.yd[n4+1]
        yv5 = domain.yd[n5+1]
        yv6 = domain.yd[n6+1]
        yv7 = domain.yd[n7+1]

        zv0 = domain.zd[n0+1]
        zv1 = domain.zd[n1+1]
        zv2 = domain.zd[n2+1]
        zv3 = domain.zd[n3+1]
        zv4 = domain.zd[n4+1]
        zv5 = domain.zd[n5+1]
        zv6 = domain.zd[n6+1]
        zv7 = domain.zd[n7+1]

        vol = domain.volo[i]*domain.vnew[i]
        norm = 1.0 / ( vol + ptiny )

        function sum4(x1,x2,x3,x4)
            x1+x2+x3+x4
        end

        dxj = -0.25*(sum4(x0,x1,x5,x4) - sum4(x3,x2,x6,x7))
        dyj = -0.25*(sum4(y0,y1,y5,y4) - sum4(y3,y2,y6,y7))
        dzj = -0.25*(sum4(z0,z1,z5,z4) - sum4(z3,z2,z6,z7))

        dxi =  0.25*(sum4(x1,x2,x6,x5) - sum4(x0,x3,x7,x4))
        dyi =  0.25*(sum4(y1,y2,y6,y5) - sum4(y0,y3,y7,y4))
        dzi =  0.25*(sum4(z1,z2,z6,z5) - sum4(z0,z3,z7,z4))

        dxk =  0.25*(sum4(x4,x5,x6,x7) - sum4(x0,x1,x2,x3))
        dyk =  0.25*(sum4(y4,y5,y6,y7) - sum4(y0,y1,y2,y3))
        dzk =  0.25*(sum4(z4,z5,z6,z7) - sum4(z0,z1,z2,z3))

        # find delvk and delxk ( i cross j )

        ax = dyi*dzj - dzi*dyj
        ay = dzi*dxj - dxi*dzj
        az = dxi*dyj - dyi*dxj

        domain.delx_zeta[i] = vol / sqrt(ax*ax + ay*ay + az*az + ptiny)

        ax = ax * norm
        ay = ay * norm
        az = az * norm

        dxv = 0.25*(sum4(xv4,xv5,xv6,xv7) - sum4(xv0,xv1,xv2,xv3))
        dyv = 0.25*(sum4(yv4,yv5,yv6,yv7) - sum4(yv0,yv1,yv2,yv3))
        dzv = 0.25*(sum4(zv4,zv5,zv6,zv7) - sum4(zv0,zv1,zv2,zv3))

        domain.delv_zeta[i] = ax*dxv + ay*dyv + az*dzv

        # find delxi and delvi ( j cross k )

        ax = dyj*dzk - dzj*dyk
        ay = dzj*dxk - dxj*dzk
        az = dxj*dyk - dyj*dxk

        domain.delx_xi[i] = vol / sqrt(ax*ax + ay*ay + az*az + ptiny)

        ax = ax * norm
        ay = ay * norm
        az = az * norm

        dxv = 0.25*(sum4(xv1,xv2,xv6,xv5) - sum4(xv0,xv3,xv7,xv4))
        dyv = 0.25*(sum4(yv1,yv2,yv6,yv5) - sum4(yv0,yv3,yv7,yv4))
        dzv = 0.25*(sum4(zv1,zv2,zv6,zv5) - sum4(zv0,zv3,zv7,zv4))

        domain.delv_xi[i] = ax*dxv + ay*dyv + az*dzv ;

        # find delxj and delvj ( k cross i )

        ax = dyk*dzi - dzk*dyi ;
        ay = dzk*dxi - dxk*dzi ;
        az = dxk*dyi - dyk*dxi ;

        domain.delx_eta[i] = vol / sqrt(ax*ax + ay*ay + az*az + ptiny)

        ax = ax * norm
        ay = ay * norm
        az = az * norm

        dxv = -0.25*(sum4(xv0,xv1,xv5,xv4) - sum4(xv3,xv2,xv6,xv7))
        dyv = -0.25*(sum4(yv0,yv1,yv5,yv4) - sum4(yv3,yv2,yv6,yv7))
        dzv = -0.25*(sum4(zv0,zv1,zv5,zv4) - sum4(zv3,zv2,zv6,zv7))

        domain.delv_eta[i] = ax*dxv + ay*dyv + az*dzv ;
    end
    return nothing
end

function calcMonotonicQRegionForElems(domain::Domain, qlc_monoq, qqc_monoq,
                                        monoq_limiter_mult, monoq_max_slope,
                                        ptiny,
                                        elength
                                    )
    for ielem in 1:elength
        i = domain.matElemlist[ielem]
        bcMask = domain.elemBC[i]

    #   phixi
        norm = 1.0 / ( domain.delv_xi[i] + ptiny )

        case = bcMask & XI_M
        if case == 0 || case == XI_M_COMM
            # MAYBE
            delvm = domain.delv_xi[domain.lxim[i]+1]
        elseif case == XI_M_SYMM
            delvm = domain.delv_xi[i]
        elseif case == XI_M_FREE
            delvm = 0.0
        else
            error("Error")
            delvm = 0.0
        end

        case = bcMask & XI_P
        if case == 0 || case == XI_P_COMM
            delvp = domain.delv_xi[domain.lxip[i]]
        elseif case == XI_P_SYMM
            delvp = domain.delv_xi[i]
        elseif case == XI_P_FREE
            delvp = 0.0
        else
            error("Error")
            delvp = 0.0
        end

        delvm = delvm * norm
        delvp = delvp * norm

        phixi = 0.5 * ( delvm + delvp )

        delvm = delvm * monoq_limiter_mult
        delvp = delvp * monoq_limiter_mult

        if  delvm < phixi
            phixi = delvm
        end
        if  delvp < phixi
            phixi = delvp
        end
        if  phixi < 0.0
            phixi = 0.0
        end
        if  phixi > monoq_max_slope
            phixi = monoq_max_slope
        end


    #   phieta
        norm = 1.0 / ( domain.delv_eta[i] + ptiny )

        case = bcMask & ETA_M
        if case == 0 || case == ETA_M_COMM
            delvm = domain.delv_eta[domain.letam[i]+1]
        elseif case == ETA_M_SYMM
            delvm = domain.delv_eta[i]
        elseif case == ETA_M_FREE
            delvm = 0.0
        else
            delvm = 0.0
            error("Error")
        end

        case = bcMask & ETA_P
        if case == 0 || case == ETA_P_COMM
            delvp = domain.delv_eta[domain.letap[i]+1]
        elseif case == ETA_P_SYMM
            delvp = domain.delv_eta[i]
        elseif case == ETA_P_FREE
            delvp = 0.0
        else
            delvp = 0.0
            error("Error")
        end

        delvm = delvm * norm
        delvp = delvp * norm

        phieta = 0.5 * ( delvm + delvp )

        delvm = delvm * monoq_limiter_mult
        delvp = delvp * monoq_limiter_mult

        if delvm  < phieta
            phieta = delvm
        end
        if delvp  < phieta
            phieta = delvp
        end
        if phieta < 0.0
            phieta = 0.0
        end
        if phieta > monoq_max_slope
            phieta = monoq_max_slope
        end

    #   phizeta
        norm = 1.0 / ( domain.delv_zeta[i] + ptiny )

        case = bcMask & ZETA_M
        if case == 0 || case == ZETA_M_COMM
            delvm = domain.delv_zeta[domain.lzetam[i]+1]
        elseif case == ZETA_M_SYMM
            delvm = domain.delv_zeta[i]
        elseif case == ZETA_M_FREE
            delvm = 0.0
        else
            delvm = 0.0
            error("Error")
        end

        case = bcMask & ZETA_P
        if case == 0 || case == ZETA_P_COMM
            delvp = domain.delv_zeta[domain.lzetap[i]+1]
        elseif case == ZETA_P_SYMM
            delvp = domain.delv_zeta[i]
        elseif case == ZETA_P_FREE
            delvp = 0.0
        else
            delvp = 0.0
            error("Error")
        end

        delvm = delvm * norm
        delvp = delvp * norm

        phizeta = 0.5 * ( delvm + delvp )

        delvm = delvm * monoq_limiter_mult
        delvp = delvp * monoq_limiter_mult

        if delvm < phizeta
            phizeta = delvm
        end
        if delvp < phizeta
            phizeta = delvp
        end
        if phizeta < 0.0
            phizeta = 0.0
        end
        if phizeta > monoq_max_slope
            phizeta = monoq_max_slope
        end

    #   Remove length scale

        if domain.vdov[i] > 0.0
            qlin  = 0.0
            qquad = 0.0
        else
            delvxxi   = domain.delv_xi[i]   * domain.delx_xi[i]
            delvxeta  = domain.delv_eta[i]  * domain.delx_eta[i]
            delvxzeta = domain.delv_zeta[i] * domain.delx_zeta[i]

            if delvxxi   > 0.0
                delvxxi   = 0.0
            end
            if delvxeta  > 0.0
                delvxeta  = 0.0
            end
            if delvxzeta > 0.0
                delvxzeta = 0.0
            end

            rho = domain.elemMass[i] / (domain.volo[i] * domain.vnew[i])

            qlin = -qlc_monoq * rho *
                    (  delvxxi   * (1.0 - phixi)  +
                        delvxeta  * (1.0 - phieta) +
                        delvxzeta * (1.0 - phizeta)  )

            qquad = qqc_monoq * rho *
                    (  delvxxi*delvxxi     * (1.0 - phixi*phixi)   +
                        delvxeta*delvxeta   * (1.0 - phieta*phieta) +
                        delvxzeta*delvxzeta * (1.0 - phizeta*phizeta)  )
        end

        domain.qq[i] = qquad
        domain.ql[i] = qlin
    end
end


function calcMonotonicQForElems(domain::Domain)

    ptiny = 1e-36
    #
    # initialize parameters
    #
    monoq_max_slope    = domain.monoq_max_slope
    monoq_limiter_mult = domain.monoq_limiter_mult

    #
    # calculate the monotonic q for pure regions
    #
    elength = domain.numElem
    if elength > 0
        qlc_monoq = domain.qlc_monoq
        qqc_monoq = domain.qqc_monoq
        calcMonotonicQRegionForElems(domain, qlc_monoq, qqc_monoq,
                                        monoq_limiter_mult,
                                        monoq_max_slope,
                                        ptiny, elength )
    end
end



function calcQForElems(domain::Domain)


    qstop = domain.qstop
    numElem = domain.numElem
    numElem8 = numElem * 8

    # MONOTONIC Q option
    commRecv(domain, MSG_MONOQ, 3,
             domain.sizeX, domain.sizeY, domain.sizeZ,
             true, true)



    # JACC Parallel


    # JACC Normal Variables
    ptiny = 1.e-36
    x_single = Vector{Float64}(undef, numElem8)
    y_single = Vector{Float64}(undef, numElem8)
    z_single = Vector{Float64}(undef, numElem8)
    xv = Vector{Float64}(undef, numElem8)
    yv = Vector{Float64}(undef, numElem8)
    zv = Vector{Float64}(undef, numElem8)

    zero_array = zeros(Float64, 10)

    # JACC Arrays
    nodelist_d = JACC.Array(domain.nodelist)
    x_d = JACC.Array(domain.x)
    y_d = JACC.Array(domain.y)
    z_d = JACC.Array(domain.z)
    xd_d = JACC.Array(domain.xd)
    yd_d = JACC.Array(domain.yd)
    zd_d = JACC.Array(domain.zd)
    volo_d = JACC.Array(domain.volo)
    vnew_d = JACC.Array(domain.vnew)
    delx_zeta_d = JACC.Array(domain.delx_zeta)
    delv_zeta_d = JACC.Array(domain.delv_zeta)
    delx_xi_d = JACC.Array(domain.delx_xi)
    delv_xi_d = JACC.Array(domain.delv_xi)
    delx_eta_d = JACC.Array(domain.delx_eta)
    delv_eta_d = JACC.Array(domain.delv_eta)

    zero_array_d = JACC.Array(zero_array)


    # JACC Big Arrays
    x_single_d = JACC.Array(x_single)
    y_single_d = JACC.Array(y_single)
    z_single_d = JACC.Array(z_single)
    xv_d = JACC.Array(xv)
    yv_d = JACC.Array(yv)
    zv_d = JACC.Array(zv)


    # time1 = time()
    # # Calculate velocity gradients
    # calcMonotonicQGradientsForElems(domain)
    # time2 = time()
    # println("calcMonotonicQGradientsForElems: ", time2-time1)

    time1 = time()
    # Calculate velocity gradients with JACC
    JACC.parallel_for(domain.numElem, MultiDcalcMonotonicQGradientsForElems, nodelist_d, x_single_d, y_single_d, z_single_d, x_d, y_d, z_d, xv_d, yv_d, zv_d, xd_d, yd_d, zd_d, volo_d, vnew_d, delx_zeta_d, delv_zeta_d, delx_xi_d, delv_xi_d, delx_eta_d, delv_eta_d, ptiny, zero_array_d)
    # JACC.parallel_for(domain.numElem, MultiDcalcMonotonicQGradientsForElems, nodelist_d, x_d, y_d, z_d, xd_d, yd_d, zd_d, volo_d, vnew_d, delx_zeta_d, delv_zeta_d, delx_xi_d, delv_xi_d, delx_eta_d, delv_eta_d, ptiny)
    time2 = time()
    println("MultiDcalcMonotonicQGradientsForElems: ", time2-time1)


    total_time = 0.0
    num_runs = 10

    for i in 1:num_runs
        time1 = time()
        # Calculate velocity gradients with JACC
        JACC.parallel_for(domain.numElem, MultiDcalcMonotonicQGradientsForElems, nodelist_d, x_single_d, y_single_d, z_single_d, x_d, y_d, z_d, xv_d, yv_d, zv_d, xd_d, yd_d, zd_d, volo_d, vnew_d, delx_zeta_d, delv_zeta_d, delx_xi_d, delv_xi_d, delx_eta_d, delv_eta_d, ptiny, zero_array_d)
        # JACC.parallel_for(domain.numElem, MultiDcalcMonotonicQGradientsForElems, nodelist_d, x_d, y_d, z_d, xd_d, yd_d, zd_d, volo_d, vnew_d, delx_zeta_d, delv_zeta_d, delx_xi_d, delv_xi_d, delx_eta_d, delv_eta_d, ptiny)
        time2 = time()
        total_time += (time2 - time1)
    end

    average_time = total_time / num_runs
    println("Average time for MultiDcalcMonotonicQGradientsForElems over $num_runs runs: ", average_time)
    # println("I'm done")
    # JACC MOVING

    # domain.delv_xi = Array(delv_xi_d)
    # domain.delv_eta = Array(delv_eta_d)
    # domain.delv_zeta = Array(delv_zeta_d)

    copyto!(domain.delx_zeta, delx_zeta_d)
    copyto!(domain.delx_xi, delx_xi_d)
    copyto!(domain.delx_eta, delx_eta_d)
    copyto!(domain.delv_xi, delv_xi_d)
    copyto!(domain.delv_eta, delv_eta_d)
    copyto!(domain.delv_zeta, delv_zeta_d)



    # Calculate velocity gradients
    # calcMonotonicQGradientsForElems(domain)

    # Transfer veloctiy gradients in the first order elements
    # problem->commElements->Transfer(CommElements::monoQ)

    fields = (domain.delv_xi, domain.delv_eta, domain.delv_zeta)
    commSend(domain, MSG_MONOQ, fields,
             domain.sizeX, domain.sizeY, domain.sizeZ,
             true, true)
    commMonoQ(domain)

    calcMonotonicQForElems(domain)

    # Don't allow excessive artificial viscosity
    if numElem != 0
        idx = -1
        for i in 1:numElem
            if domain.q[i] > qstop
                idx = i
                break
            end
        end

        if idx >= 0
            if domain.comm !== nothing
                MPI.Abort(MPI.COMM_WORLD, 1)
            end
            error("QStopError")
        end
    end
end


@inline function MultiDcalcPressureForElems(i, p_new, bvc, pbvc, e_old, compression, vnewc, pmin, p_cut, eosvmax)
    @inbounds begin
        c1s = 2.0 / 3.0

        # Calculate bvc and pbvc
        bvc[i] = c1s * (compression[i] + 1.0)
        pbvc[i] = c1s

        # Calculate new pressure
        p_new[i] = bvc[i] * e_old[i]

        # Apply pressure cut-off
        if abs(p_new[i]) < p_cut
            p_new[i] = 0.0
        end

        # Apply volume condition
        if vnewc[i] >= eosvmax
            p_new[i] = 0.0
        end

        # Apply minimum pressure condition
        if p_new[i] < pmin
            p_new[i] = pmin
        end
    end
end


function calcPressureForElems(domain::Domain, p_new, bvc,
                                 pbvc, e_old,
                                 compression, vnewc,
                                 pmin,
                                 p_cut,eosvmax,
                                 length              )

    c1s = 2.0/3.0

    for i in 1:length
        bvc[i] = c1s * (compression[i] + 1.0)
        pbvc[i] = c1s
    end

    for i in 1:length
        p_new[i] = bvc[i] * e_old[i]

        if abs(p_new[i]) < p_cut
            p_new[i] = 0.0
        end

        if vnewc[i] >= eosvmax # impossible condition here?
            p_new[i] = 0.0
        end

        if p_new[i] < pmin
            p_new[i] = pmin
        end
    end
end


function MultiDcalcEnergyForElems(i, p_new,  e_new,  q_new,
                                bvc,  pbvc,
                                p_old,  e_old,  q_old,
                                compression,  compHalfStep,
                                vnewc,  work,  delvc,  pmin,
                                p_cut,   e_cut,  q_cut,  emin,
                                qq,  ql,
                                rho0,
                                eosvmax,
                                pHalfStep                     )

    TINY1 = 0.111111e-36
    TINY3 = 0.333333e-18
    SIXTH = 1.0 / 6.0

    if(i > 27000)
        return
    end

    e_new[i] = e_old[i] - 0.5 * delvc[i] * (p_old[i] + q_old[i]) + 0.5 * work[i]

    if e_new[i]  < emin
        e_new[i] = emin
    end

    MultiDcalcPressureForElems(i, pHalfStep, bvc, pbvc, e_new, compHalfStep,
                            vnewc, pmin, p_cut, eosvmax)

    vhalf = 1.0 / (1.0 + compHalfStep[i])

    if  delvc[i] > 0.0
        q_new[i] = 0.0
    else
        ssc = (( pbvc[i] * e_new[i]
            + vhalf * vhalf * bvc[i] * pHalfStep[i] ) / rho0)
        if ssc <= TINY1
            ssc = TINY3
        else
            ssc = sqrt(ssc)
        end

        q_new[i] = (ssc*ql[i] + qq[i])
    end

    e_new[i] = (e_new[i] + 0.5 * delvc[i] * (  3.0*(p_old[i]     + q_old[i])
        - 4.0*(pHalfStep[i] + q_new[i])))
    
    e_new[i] = e_new[i] + 0.5 * work[i]
    if abs(e_new[i]) < e_cut
        e_new[i] = 0.0
    end
    if e_new[i]  < emin
        e_new[i] = emin
    end

    MultiDcalcPressureForElems(i, p_new, bvc, pbvc, e_new, compression,
                            vnewc, pmin, p_cut, eosvmax)

    # for i in 1:length
    if delvc[i] > 0.0
        q_tilde = 0.0
    else
        ssc = ( pbvc[i] * e_new[i]
            + vnewc[i] * vnewc[i] * bvc[i] * p_new[i] ) / rho0

        if ssc <= TINY1
            ssc = TINY3
        else
            ssc = sqrt(ssc)
        end

        q_tilde = (ssc*ql[i] + qq[i])
    end

    e_new[i] = (e_new[i] - (  7.0*(p_old[i]     + q_old[i])
                        -    8.0*(pHalfStep[i] + q_new[i])
                        + (p_new[i] + q_tilde)) * delvc[i]*SIXTH)

    if abs(e_new[i]) < e_cut
        e_new[i] = 0.0
    end
    if e_new[i]  < emin
        e_new[i] = emin
    end
    # end

    MultiDcalcPressureForElems(i, p_new, bvc, pbvc, e_new, compression,
                            vnewc, pmin, p_cut, eosvmax)

    # for i in 1:length

    if delvc[i] <= 0.0
        ssc = (( pbvc[i] * e_new[i]
            + vnewc[i] * vnewc[i] * bvc[i] * p_new[i] ) / rho0)

        if ssc <= TINY1
            ssc = TINY3
        else
            ssc = sqrt(ssc)
        end

        q_new[i] = (ssc*ql[i] + qq[i])

        if abs(q_new[i]) < q_cut
            q_new[i] = 0.0
        end
    end
    # end
end


function calcEnergyForElems(domain::Domain, p_new,  e_new,  q_new,
                                bvc,  pbvc,
                                p_old,  e_old,  q_old,
                                compression,  compHalfStep,
                                vnewc,  work,  delvc,  pmin,
                                p_cut,   e_cut,  q_cut,  emin,
                                qq,  ql,
                                rho0,
                                eosvmax,
                                length                          )

    TINY1 = 0.111111e-36
    TINY3 = 0.333333e-18
    SIXTH = 1.0 / 6.0


    pHalfStep = Vector{Float64}(undef, length)

    for i in 1:length
        e_new[i] = e_old[i] - 0.5 * delvc[i] * (p_old[i] + q_old[i]) + 0.5 * work[i]

        if e_new[i]  < emin
            e_new[i] = emin
        end
    end

    calcPressureForElems(domain, pHalfStep, bvc, pbvc, e_new, compHalfStep,
                            vnewc, pmin, p_cut, eosvmax, length)
    for i in 1:length
        vhalf = 1.0 / (1.0 + compHalfStep[i])

        if  delvc[i] > 0.0
        #      q_new(i) /* = qq(i) = ql(i) */ = Real_t(0.) ;
            q_new[i] = 0.0
        else
            ssc = (( pbvc[i] * e_new[i]
                + vhalf * vhalf * bvc[i] * pHalfStep[i] ) / rho0)
            if ssc <= TINY1
                ssc = TINY3
            else
                ssc = sqrt(ssc)
            end

            q_new[i] = (ssc*ql[i] + qq[i])
        end

        e_new[i] = (e_new[i] + 0.5 * delvc[i] * (  3.0*(p_old[i]     + q_old[i])
            - 4.0*(pHalfStep[i] + q_new[i])))
    end

    for i in 1:length
        e_new[i] = e_new[i] + 0.5 * work[i]
        if abs(e_new[i]) < e_cut
            e_new[i] = 0.0
        end
        if e_new[i]  < emin
            e_new[i] = emin
        end
    end

    calcPressureForElems(domain, p_new, bvc, pbvc, e_new, compression,
                            vnewc, pmin, p_cut, eosvmax, length)

    for i in 1:length
        if delvc[i] > 0.0
            q_tilde = 0.0
        else
            ssc = ( pbvc[i] * e_new[i]
                + vnewc[i] * vnewc[i] * bvc[i] * p_new[i] ) / rho0

            if ssc <= TINY1
                ssc = TINY3
            else
                ssc = sqrt(ssc)
            end

            q_tilde = (ssc*ql[i] + qq[i])
        end

        e_new[i] = (e_new[i] - (  7.0*(p_old[i]     + q_old[i])
                            -    8.0*(pHalfStep[i] + q_new[i])
                            + (p_new[i] + q_tilde)) * delvc[i]*SIXTH)

        if abs(e_new[i]) < e_cut
            e_new[i] = 0.0
        end
        if e_new[i]  < emin
            e_new[i] = emin
        end
    end

    calcPressureForElems(domain, p_new, bvc, pbvc, e_new, compression,
                            vnewc, pmin, p_cut, eosvmax, length)

    for i in 1:length

        if delvc[i] <= 0.0
            ssc = (( pbvc[i] * e_new[i]
                + vnewc[i] * vnewc[i] * bvc[i] * p_new[i] ) / rho0)

            if ssc <= TINY1
                ssc = TINY3
            else
                ssc = sqrt(ssc)
            end

            q_new[i] = (ssc*ql[i] + qq[i])

            if abs(q_new[i]) < q_cut
                q_new[i] = 0.0
            end
        end
    end
end

function calcSoundSpeedForElems(domain::Domain, vnewc,  rho0, enewc,
                                  pnewc, pbvc,
                                  bvc, ss4o3, nz       )
    TINY1 = 0.111111e-36
    TINY3 = 0.333333e-18

    for i in 1:nz
        iz = domain.matElemlist[i]
        ssTmp = ((pbvc[i] * enewc[i] + vnewc[i] * vnewc[i] *
                            bvc[i] * pnewc[i]) / rho0)
        if ssTmp <= TINY1
            ssTmp = TINY3
        else
            ssTmp = sqrt(ssTmp)
        end
        domain.ss[iz] = ssTmp
    end
end


function evalEOSForElems(domain::Domain, vnewc, length)

    e_cut = domain.e_cut
    p_cut = domain.p_cut
    ss4o3 = domain.ss4o3
    q_cut = domain.q_cut

    eosvmax = domain.eosvmax
    eosvmin = domain.eosvmin
    pmin    = domain.pmin
    emin    = domain.emin
    rho0    = domain.refdens

    e_old = Vector{Float64}(undef, length)
    delvc = Vector{Float64}(undef, length)
    p_old = Vector{Float64}(undef, length)
    q_old = Vector{Float64}(undef, length)
    compression = Vector{Float64}(undef, length)
    compHalfStep = Vector{Float64}(undef, length)
    qq = Vector{Float64}(undef, length)
    ql = Vector{Float64}(undef, length)
    work = Vector{Float64}(undef, length)
    p_new = Vector{Float64}(undef, length)
    e_new = Vector{Float64}(undef, length)
    q_new = Vector{Float64}(undef, length)
    bvc = Vector{Float64}(undef, length)
    pbvc = Vector{Float64}(undef, length)

    # compress data, minimal set
    for i in 1:length
        zidx = domain.matElemlist[i]
        e_old[i] = domain.e[zidx]
    end

    for i in 1:length
        zidx = domain.matElemlist[i]
        delvc[i] = domain.delv[zidx]
    end

    for i in 1:length
        zidx = domain.matElemlist[i]
        p_old[i] = domain.p[zidx]
    end

    for i in 1:length
        zidx = domain.matElemlist[i]
        q_old[i] = domain.q[zidx]
    end

    for i in 1:length
        compression[i] = 1.0 / vnewc[i] - 1.0
        vchalf = vnewc[i] - delvc[i] * 0.5
        compHalfStep[i] = (1.0 / vchalf) - 1.0
    end

    # Check for v > eosvmax or v < eosvmin
    if eosvmin != 0.0
        for i in 1:length
            if vnewc[i] <= eosvmin  # impossible due to calling func?
                compHalfStep[i] = compression[i]
            end
        end
    end

    if eosvmax != 0.0
        for i in 1:length
            if vnewc[i] >= eosvmax  # impossible due to calling func?
                p_old[i]        = 0.0
                compression[i]  = 0.0
                compHalfStep[i] = 0.0
            end
        end
    end
    for i in 1:length
        zidx = domain.matElemlist[i]
        qq[i] = domain.qq[zidx]
        ql[i] = domain.ql[zidx]
        work[i] = 0.0
    end


    pHalfStep = Vector{Float64}(undef, length)

    p_new_d = JACC.Array(p_new)
    e_new_d = JACC.Array(e_new)
    q_new_d = JACC.Array(q_new)
    bvc_d = JACC.Array(bvc)
    pbvc_d = JACC.Array(pbvc)
    p_old_d = JACC.Array(p_old)
    e_old_d = JACC.Array(e_old)
    q_old_d = JACC.Array(q_old)
    compression_d = JACC.Array(compression)
    compHalfStep_d = JACC.Array(compHalfStep)
    vnewc_d = JACC.Array(vnewc)
    work_d = JACC.Array(work)
    delvc_d = JACC.Array(delvc)
    qq_d = JACC.Array(qq)
    ql_d = JACC.Array(ql)
    pHalfStep_d = JACC.Array(pHalfStep)

    time1 = time()

    JACC.parallel_for(length, MultiDcalcEnergyForElems, p_new_d, e_new_d, q_new_d,
                        bvc_d, pbvc_d, p_old_d, e_old_d, q_old_d, compression_d,
                        compHalfStep_d, vnewc_d, work_d, delvc_d, pmin,
                        p_cut, e_cut, q_cut, emin, qq_d, ql_d, rho0, eosvmax, pHalfStep_d)
    # CUDA.synchronize()
    time2 = time()

    println("Time for MultiDcalcEnergyForElems: ", time2 - time1)

    total_time = 0.0
    num_runs = 10

    for jk in 1:num_runs
        time1 = time()

        JACC.parallel_for(length, MultiDcalcEnergyForElems, p_new_d, e_new_d, q_new_d,
                            bvc_d, pbvc_d, p_old_d, e_old_d, q_old_d, compression_d,
                            compHalfStep_d, vnewc_d, work_d, delvc_d, pmin,
                            p_cut, e_cut, q_cut, emin, qq_d, ql_d, rho0, eosvmax, pHalfStep_d)
        # CUDA.synchronize()
        time2 = time()

        total_time = total_time + (time2 - time1)
    end

    average_time = total_time / num_runs
    println("Average time for MultiDcalcEnergyForElems over $num_runs runs: ", average_time)

    # Use these with JACC

    p_new = Array(p_new_d)
    e_new = Array(e_new_d)
    q_new = Array(q_new_d)
    bvc = Array(bvc_d)
    pbvc = Array(pbvc_d)
    vnewc = Array(vnewc_d)


    for i in 1:length
        zidx = domain.matElemlist[i]
        domain.p[zidx] = p_new[i]
    end

    for i in 1:length
        zidx = domain.matElemlist[i]
        domain.e[zidx] = e_new[i]
    end

    for i in 1:length
        zidx = domain.matElemlist[i]
        domain.q[zidx] = q_new[i]
    end


    calcSoundSpeedForElems(domain, vnewc, rho0, e_new, p_new,
                                pbvc, bvc, ss4o3, length)
end


function applyMaterialPropertiesForElems(domain::Domain)

    length = domain.numElem

    if length != 0
        #  Expose all of the variables needed for material evaluation
        eosvmin = domain.eosvmin
        eosvmax = domain.eosvmax
        vnewc = Vector{Float64}(undef, length)

        for i in 1:length
            zn = domain.matElemlist[i]
            vnewc[i] = domain.vnew[zn]
        end

        if eosvmin != 0.0
            for i in 1:length
                if vnewc[i] < eosvmin
                    vnewc[i] = eosvmin
                end
            end
        end

        if eosvmax != 0.0
            for i in 1:length
                if vnewc[i] > eosvmax
                    vnewc[i] = eosvmax
                end
            end
        end

        for i in 1:length
            zn = domain.matElemlist[i]
            vc = domain.v[zn]
            if eosvmin != 0.0
                if vc < eosvmin
                    vc = eosvmin
                end
            end
            if eosvmax != 0.0
                if vc > eosvmax
                    vc = eosvmax
                end
            end
            if vc <= 0.0
                error("Volume Error :2887")
            end
        end
        evalEOSForElems(domain::Domain, vnewc, length)
    end
end

function updateVolumesForElems(domain::Domain)
    numElem = domain.numElem

    if numElem != 0
        v_cut = domain.v_cut

        for i in 1:numElem
            tmpV = domain.vnew[i]

            if abs(tmpV - 1.0) < v_cut
                tmpV = 1.0
            end
            domain.v[i] = tmpV
            if tmpV <= 0.0
                error("Volume Error :2908")
            end
        end
    end
end

function lagrangeElements(domain::Domain)

    delt = domain.deltatime
    domain.vnew = Vector{Float64}(undef, domain.numElem)
    domain.dxx = Vector{Float64}(undef, domain.numElem)
    domain.dyy = Vector{Float64}(undef, domain.numElem)
    domain.dzz = Vector{Float64}(undef, domain.numElem)

    domain.delx_xi = Vector{Float64}(undef, domain.numElem)
    domain.delx_eta = Vector{Float64}(undef, domain.numElem)
    domain.delx_zeta = Vector{Float64}(undef, domain.numElem)

    allElem = domain.numElem +  # local elem
            2*domain.sizeX*domain.sizeY + # plane ghosts
            2*domain.sizeX*domain.sizeZ + # row ghosts
            2*domain.sizeY*domain.sizeZ

    domain.delv_xi = Vector{Float64}(undef, allElem)
    domain.delv_eta = Vector{Float64}(undef, allElem)
    domain.delv_zeta = Vector{Float64}(undef, allElem)

    calcLagrangeElements(domain, delt)

    # Calculate Q.  (Monotonic q option requires communication)
    calcQForElems(domain)

    applyMaterialPropertiesForElems(domain)

    updateVolumesForElems(domain)

end

function calcCourantConstraintForElems(domain::Domain)

    dtcourant    = 1.0e+20
    courant_elem = -1

    qqc = domain.qqc
    length = domain.numElem

    qqc2 = 64.0 * qqc * qqc

    # Rewritten OpenMP code to sequential code
    courant_elem_per_thread = -1
    dtcourant_per_thread =  1.0e+20


    for i in 1:length
        indx = domain.matElemlist[i]

        dtf = domain.ss[indx] * domain.ss[indx]

        if domain.vdov[indx] < 0.0

        dtf = (dtf + qqc2 * domain.arealg[indx] * domain.arealg[indx]
                    * domain.vdov[indx]* domain.vdov[indx])
        end

        dtf = sqrt(dtf)

        dtf = domain.arealg[indx] / dtf

        #  determine minimum timestep with its corresponding elem
        if domain.vdov[indx] != 0.0
            if dtf < dtcourant_per_thread

                dtcourant_per_thread = dtf
                courant_elem_per_thread = indx
            end
        end
    end

    if dtcourant_per_thread < dtcourant
        dtcourant = dtcourant_per_thread
        courant_elem =  courant_elem_per_thread
    end


    # Don't try to register a time constraint if none of the elements
    # were active
    if courant_elem != -1
        domain.dtcourant = dtcourant
    end

    return nothing
end


function calcHydroConstraintForElems(domain::Domain)

    dthydro = 1.0e+20
    hydro_elem = -1
    dvovmax = domain.dvovmax
    length = domain.numElem

    # Rewritten OpenMP code to sequential code

    hydro_elem_per_thread = hydro_elem
    dthydro_per_thread = dthydro

    for i in 1:length
        indx = domain.matElemlist[i]

        if domain.vdov[indx] != 0.0
            dtdvov = dvovmax / (abs(domain.vdov[indx])+1.e-20)

            if dthydro_per_thread > dtdvov
                dthydro_per_thread = dtdvov
                hydro_elem_per_thread = indx
            end
        end
    end

    if dthydro_per_thread < dthydro
      dthydro = dthydro_per_thread
      hydro_elem =  hydro_elem_per_thread
    end

    if hydro_elem != -1
        domain.dthydro = dthydro
    end
    return nothing
end



function calcTimeConstraintsForElems(domain::Domain)
  # evaluate time constraint
  calcCourantConstraintForElems(domain::Domain)

  # check hydro constraint
  calcHydroConstraintForElems(domain::Domain)
end


function lagrangeLeapFrog(domain::Domain)

   # calculate nodal forces, accelerations, velocities, positions, with
   # applied boundary conditions and slide surface considerations */
   # Time increment
   lagrangeNodal(domain)

   # calculate element quantities (i.e. velocity gradient & q), and update
   # material states */
   lagrangeElements(domain)

   calcTimeConstraintsForElems(domain)
   return nothing
end
