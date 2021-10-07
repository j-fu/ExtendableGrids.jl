
struct CellFinder{Tv,Ti}
    xgrid::ExtendableGrid{Tv,Ti}
    xCellFaces::Adjacency{Ti}
    xFaceCells::Adjacency{Ti}
    xCellGeometries::GridEGTypes
    facetogo::Array{Array{Ti,1},1}
    previous_cells::Array{Ti,1}
    EG::Array{Type{<:AbstractElementGeometry},1}
    L2G4EG::Vector{L2GTransformer{Tv,Ti,EG,CS} where {EG<:AbstractElementGeometry,CS<:AbstractCoordinateSystem}} # each Geometry has its own L2GTransformer
    invA::Matrix{Tv}
    xreftest::Array{Tv,1}
    cx::Vector{Tv}
end

## postprocessing function that calcualte conditions to decide over which face to go next
function postprocess_xreftest!(xreftest::Array{Tv}, ::Type{<:Edge1D}) where {Tv}
    xreftest[2] = 1 - xreftest[1]
    return nothing
end
function postprocess_xreftest!(xreftest::Array{Tv}, ::Type{<:Triangle2D}) where {Tv}
    xreftest[3] = 1 - xreftest[1] - xreftest[2]
    return nothing
end
function postprocess_xreftest!(xreftest::Array{Tv}, ::Type{<:Tetrahedron3D}) where {Tv}
    xreftest[4] = 1 - xreftest[1] - xreftest[2] - xreftest[3]
    return nothing
end
function postprocess_xreftest!(xreftest::Array{Tv}, ::Type{<:Parallelogram2D}) where {Tv}
    xreftest[3] = 1 - xreftest[1]
    xreftest[4] = 1 - xreftest[2]
    return nothing
end
function postprocess_xreftest!(xreftest::Array{Tv}, ::Type{<:Parallelepiped3D}) where {Tv}
    xreftest[4] = 1 - xreftest[1]
    xreftest[5] = 1 - xreftest[2]
    xreftest[6] = 1 - xreftest[3]
    return nothing
end


function CellFinder(xgrid::ExtendableGrid{Tv,Ti}) where {Tv,Ti}
    CS = xgrid[CoordinateSystem]
    EG = xgrid[UniqueCellGeometries]
    L2G4EG = Vector{L2GTransformer{Tv,Ti,EG,CS} where {EG<:AbstractElementGeometry,CS<:AbstractCoordinateSystem}}(undef,length(EG))
    facetogo = Array{Array{Ti,1},1}(undef,length(EG))
    xreftestlength::Int = 0
    for j = 1 : length(EG)
        L2G4EG[j] = L2GTransformer(EG[j], xgrid, ON_CELLS)
        if EG[j] <: AbstractElementGeometry1D
            facetogo[j] = [1, 2]
            xreftestlength = max(xreftestlength,2)
        elseif EG[j] <: Triangle2D
            facetogo[j] = [3, 1, 2]
            xreftestlength = max(xreftestlength,3)
        elseif EG[j] <: Tetrahedron3D
            facetogo[j] = [4, 2, 1, 3]
            xreftestlength = max(xreftestlength,4)
        elseif EG[j] <: Parallelogram2D
            facetogo[j] = [4, 1, 2, 3]
            xreftestlength = max(xreftestlength,4)
        elseif EG[j] <: Parallelepiped3D
            facetogo[j] = [5, 2, 1, 3, 4, 6]
            xreftestlength = max(xreftestlength,6)
        else
            @error "ElementGeometry not supported by CellFinder"
        end
    end

    edim = dim_element(EG[1])
    A = zeros(Tv,edim,edim)
    xreftest = zeros(Tv,xreftestlength)
    return CellFinder{Tv,Ti}(xgrid, xgrid[CellFaces], xgrid[FaceCells], xgrid[CellGeometries], facetogo, zeros(Ti,3), EG, L2G4EG, A, xreftest, zeros(Tv,edim))
end


# modification of pdelib function in gfind.cxx
function gFindLocal!(xref, CF::CellFinder{Tv,Ti}, x; icellstart::Int = 1, eps = 1e-14) where{Tv,Ti}

    # works for convex domainsand simplices only !
    xCellFaces::Adjacency{Ti} = CF.xCellFaces
    xFaceCells::Adjacency{Ti} = CF.xFaceCells
    xCellGeometries::GridEGTypes = CF.xCellGeometries
    EG::GridEGTypes = CF.EG
    cx::Vector{Tv} = CF.cx
    cEG::Int = 0
    facetogo::Array{Array{Ti,1},1} = CF.facetogo
    icell::Int = icellstart
    previous_cells::Array{Ti,1} = CF.previous_cells
    fill!(previous_cells,0)
    xreftest::Array{Tv,1} = CF.xreftest
    L2G::L2GTransformer{Tv,Ti} = CF.L2G4EG[1]
    L2Gb::Vector{Tv} = L2G.b

    invA::Matrix{Tv} = CF.invA
    imin::Int = 0

    while (true)
        # find current cell geometry index
        cEG = 1
        while xCellGeometries[icell] != EG[cEG]
            cEG += 1
        end

        # update local 2 global map
        L2G = CF.L2G4EG[cEG]
        update_trafo!(L2G, icell)
        L2Gb = L2G.b

        # compute barycentric coordinates of node
        for j = 1 : length(x)
            cx[j] = x[j] - L2Gb[j]
        end
        mapderiv!(invA,L2G,xref)
        fill!(xreftest,0)
        for j = 1 : length(x), k = 1 : length(x)
            xreftest[k] += invA[j,k] * cx[j]
        end
        postprocess_xreftest!(xreftest,xCellGeometries[icell])

        # find minimal barycentric coordinate with
        imin = 1
        for i = 2 : length(facetogo[cEG])
            if xreftest[imin] >= xreftest[i]
                imin = i
            end
        end

        # if all barycentric coordinates are within [0,1] the including cell is found
        if xreftest[imin] >= -eps
            xref .= view(xreftest,1:length(xref))
            return icell
        end

        # otherwise: go into direction of minimal barycentric coordinates
        for j = 1 : length(previous_cells)-1
            previous_cells[j] = previous_cells[j+1]
        end
        previous_cells[end] = icell
        icell = xFaceCells[1,xCellFaces[facetogo[cEG][imin],icell]]
        if icell == previous_cells[end]
            icell = xFaceCells[2,xCellFaces[facetogo[cEG][imin],icell]]
            if icell == 0
                @debug  "could not find point in any cell and ended up at boundary of domain (maybe x lies outside of the domain ?)"
                return 0
            end
        end

        if icell == previous_cells[end-1]
            @debug  "could not find point in any cell and ended up going in circles (better try brute force search)"
            return 0
        end
    end
    
    return 0
end

function gFindBruteForce!(xref, CF::CellFinder{Tv,Ti}, x; eps = 1e-14) where {Tv,Ti}

    cx::Vector{Tv} = CF.cx
    cEG::Int = 0
    EG::GridEGTypes = CF.EG
    xreftest::Array{Tv,1} = CF.xreftest
    xCellGeometries::GridEGTypes = CF.xCellGeometries
    facetogo::Array{Array{Ti,1},1} = CF.facetogo
    L2Gb::Vector{Tv} = CF.L2G4EG[1].b
    invA::Matrix{Tv} = CF.invA
    imin::Int = 0

    for icell = 1 : num_sources(CF.xgrid[CellNodes])

        # find current cell geometry index
        cEG = 1
        while xCellGeometries[icell] != EG[cEG]
            cEG += 1
        end

        # update local 2 global map
        L2G = CF.L2G4EG[cEG]
        update_trafo!(L2G, icell)
        L2Gb = L2G.b

        # compute barycentric coordinates of node
        for j = 1 : length(x)
            cx[j] = x[j] - L2Gb[j]
        end
        mapderiv!(invA,L2G,xref)
        fill!(xreftest,0)
        for j = 1 : length(x), k = 1 : length(x)
            xreftest[k] += invA[j,k] * cx[j]
        end
        postprocess_xreftest!(xreftest,xCellGeometries[icell])

        # find minimal barycentric coordinate with
        imin = 1
        for i = 2 : length(facetogo[cEG])
            if xreftest[imin] >= xreftest[i]
                imin = i
            end
        end

        # if all barycentric coordinates are within [0,1] the including cell is found
        if xreftest[imin] >= -eps
            xref .= view(xreftest,1:length(xref))
            return icell
        end
    end

    @debug "gFindBruteForce did not find any cell that contains x = $x (make sure that x is inside the domain, or try reducing $eps)"
    
    return 0
end
