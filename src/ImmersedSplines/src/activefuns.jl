export spline_extension_operator, Perimeter, GenPerimeter

"""
    spline_extension_operator(U::TensorProduct{2, <:SplineSpace}, F)

Compute an extension operator that stabilizes the splinespace according to the
definition of [Höllig, Klaus, Ulrich Reif, and Joachim Wipper. "Weighted extended 
B-spline approximation of Dirichlet problems." SIAM Journal on Numerical Analysis 39, 
no. 2 (2001): 442-462.]
"""
function spline_extension_operator(U::TensorProduct{Dim, <:SplineSpace}, phi) where Dim

    # Construct closest index array classification
    ExtensionArray = ClosestExtensionArray(U, phi)

    # main dimensions
    n = map(s -> dimsplinespace(s), U)
    N = prod(n)

    # allocate extraction operator
    C = spzeros(Float64, N, N)

    # compute rows of extension operator
    for mj in ExtensionArray.perimeter.i

        row = ExtensionArray.perimeter.l[mj]

        # Loop over boundary splines
        if ExtensionArray.functions[mj.I...] == 1

            # compute closest extension indices
            mi = ExtensionArray(mj)

            # compute extension coefficients
            eᵢⱼ = KroneckerProduct(k -> compute_extension_coefficients(U[k], mi.indices[k][end], mj[k]), Dim:-1:1)

            cols = ExtensionArray.perimeter.l[mi]
            C[row, cols] = eᵢⱼ  
        
         # skip in-active splines
        elseif ExtensionArray.functions[mj.I...] != 0
            C[row,row] = 1.0
        end
    end

    # remove columns with only zeros
    colptr = unique(C.colptr)
    rowval = C.rowval
    nzval = C.nzval
    m = C.m
    n = length(colptr)-1

    # set numerical zeros to zero
    for i in eachindex(nzval)
        if abs(nzval[i]) < 100*eps()
            nzval[i] = 0
        end
    end

    # create sparse matrix
    C = SparseMatrixCSC{Float64,Int64}(m, n, colptr, rowval, nzval)
    dropzeros!(C)
    return C
end

function spline_extension_operator(U::SplineSpace, phi)
    return spline_extension_operator(⨷(U), phi)
end

function compute_extension_coefficients(S::SplineSpace, span, k)
    J = S.C.colptr[k]:S.C.colptr[k+1]-1
    I = S.C.rowval[J]
    V = S.C.nzval[J]
    J2 = S.C.colptr[span]:S.C.colptr[span+1]-1
    c = zeros(S.p+1)
    for l in eachindex(I)
        c += V[l] * compute_extension_coefficients(S.p, S.U, S.C.rowval[J2[end]], I[l])
    end
    return c ./ sum(V)
end

"""
    compute_extension_coefficients(p, kts, span, k)

Compute the extension coefficients using dual functionals. This works for
general non-uniform knot vectors and is based on the paper [Höllig, Klaus, 
and Ulrich Reif. "Nonuniform web-splines." Computer Aided Geometric Design 
20, no. 5 (2003): 277-294.]
"""
function compute_extension_coefficients(p::Degree, kts, span::Int, k::Int)
    t = sum(kts[k+1:k+p]) / p # greville point
    γ = UnivariateSplines.deboor_fix_dual_functional_coeffs(kts[k:k+p+1], t)
    B = dersbsplinebasisfuns(p, kts, span, t, p+1)
    return B * γ
end

"""
    GeneratePerimeter{Dim,Indices}

Datastructure that allows iteration over active neighborhoods of boundary 
functions
"""
struct GeneratePerimeter{Dim,Indices}
    p::NTuple{Dim, Degree}
    i::CartesianIndices{Dim,Indices}
    l::LinearIndices{Dim,Indices}
end

function GeneratePerimeter(U::TensorProduct{Dim,<:SplineSpace}) where Dim
    p = ntuple(k -> Degree(U[k]), Dim)
    i = CartesianIndices(ntuple(k -> Base.OneTo(dimsplinespace(U[k])), Dim))  
    l = LinearIndices(i)  
    return GeneratePerimeter(p, i, l)  
end

function (gperm::GeneratePerimeter{Dim})(mi::CartesianIndex{Dim}, ring::Int=1) where Dim
    p = gperm.p
    s = size(gperm.i)
    i = ntuple(k -> max(1,mi[k]-p[k]-ring):1:min(s[k]-p[k],mi[k]+ring), Dim)
    return IPerm(ring, p, mi, view(gperm.i, i...), gperm.i)
end


"""
    ClosestExtensionArray(F, U::TensorProduct{Dim,<:SplineSpace})

Defines a `closest` index array, according to the definition of 
[Höllig, Klaus, Ulrich Reif, and Joachim Wipper. "Weighted extended B-spline
approximation of Dirichlet problems." SIAM Journal on Numerical Analysis 39, 
no. 2 (2001): 442-462.]
"""
struct ClosestExtensionArray{Dim, Indices}
    points::Array{Bool,Dim}     # active breakpoints
    elements::Array{Int8,Dim}   # active elements
    functions::Array{Int8,Dim}  # active functions
    perimeter::GeneratePerimeter{Dim,Indices}
end

function ClosestExtensionArray(U::TensorProduct{Dim,<:SplineSpace}, phi) where Dim
    partition = CartesianProduct(u -> breakpoints(u), U)
    active_points = point_is_inside(phi, partition)
    active_elements = element_is_inside(active_points)
    active_functions = active_splines(U, active_elements)
    perimeter = GeneratePerimeter(U)
    return ClosestExtensionArray(active_points, active_elements, active_functions, perimeter)
end

"""
    ClosestExtensionArray(F, U::TensorProduct{Dim,<:SplineSpace})

Computes a `closest` index array, according to the definition of 
[Höllig, Klaus, Ulrich Reif, and Joachim Wipper. "Weighted extended B-spline
approximation of Dirichlet problems." SIAM Journal on Numerical Analysis 39, 
no. 2 (2001): 442-462.]
"""
function (extop::ClosestExtensionArray{Dim})(index::CartesianIndex{Dim}) where {Dim}
    return find_extension_indices(extop.functions, extop.perimeter, index)
end

# function (extop::ClosestExtensionArray{Dim})(index::CartesianIndex{Dim}) where {Dim}
#     ei = find_first_active_element(extop.active_elements, extop.perimeter) # first active element
# end

struct Perimeter{Dim, X<:CartesianIndices{Dim}}
    maxsize::NTuple{Dim,Int}
    inner::X
    outer::X
end

function Perimeter(inner::UnitRange{Int}...; maxsize)
    outer = map(x -> growone(x, 10), inner)
    return Perimeter(maxsize, CartesianIndices(inner), CartesianIndices(outer))
end

function perm_iter_imp(iter::Perimeter, state)
    while state<length(iter.outer)
        state+=1
        x = iter.outer[state]
        if x ∉ iter.inner
            return (x, state)
        end
    end
end

function Base.iterate(iter::Perimeter)
    return perm_iter_imp(iter, 0)
end

function Base.iterate(iter::Perimeter, state)
    return perm_iter_imp(iter, state)
end

function growone(x::UnitRange{Int}, maxsize)::UnitRange{Int}        
    @assert x[end]<=maxsize
    a = x[1]==1 ? 1 : x[1]-1
    b = x[end]==maxsize ? maxsize : x[end]+1
    return a:b
end

struct GenPerimeter{Dim, X<:CartesianProduct{Dim}}
    supports::X
end

function (perm::GenPerimeter)(idx::CartesianIndex{Dim}) where Dim
    α = ntuple(k -> perm.supports.data[k][idx[k]], Dim)
    m = ntuple(k -> perm.supports.data[k][end][end], Dim)
    return Perimeter(α...; maxsize=m)
end

"""
    IPerm{Dim,Indices}

Iterator that returns active neighborhoods around a B-spline.
"""
struct IPerm{Dim,Indices}
    k::Int                   # k-ring 
    p::NTuple{Dim, Int}      # polynomial degree
    c::CartesianIndex{Dim}   # current cartesian index
    i::Indices               # view of parent
    parent::CartesianIndices{Dim,NTuple{Dim,Base.OneTo{Int64}}}
end

function iter_imp(iter::IPerm{Dim}, state) where Dim
    while state<length(iter.i)
        state+=1
        mi = iter.i[state]
        for k in Dim:-1:1
            if mi[k]==iter.c[k]-iter.p[k]-iter.k || mi[k]==iter.c[k]+iter.k
                res = view(iter.parent, ntuple(k->mi[k]:mi[k]+iter.p[k], Dim)...)
                return (res, state)
            end
        end
    end
end

function Base.iterate(iter::IPerm{Dim}) where Dim
    return iter_imp(iter, 0)
end

function Base.iterate(iter::IPerm{Dim}, state) where Dim
    return iter_imp(iter, state)
end

"""
point_is_inside(phi, partition::CartesianProduct)

Return a boolean array that specifies whether a point is inside or 
outside, respectively
"""
function point_is_inside(phi, X)
    Y = zeros(Bool, size(X))
    for k in eachindex(X)
        Y[k] = phi(SVector(X[k])) <= 0.0
    end
    return Y
end


"""
element_is_inside(Y::Matrix{Bool})

Return a boolean array that specifies whether an element is inside or 
outside, respectively
"""
function element_is_inside(phi, X)
    Y = point_is_inside(phi, X)
    return element_is_inside(Y)
end

function element_is_inside(Y::Vector{Bool})
    m = map(s -> s-1, length(Y))
    E = zeros(Int8, m)
    for i in eachindex(E)
        E[i] = sum(Y[i:i+1])
    end
    return E
end    

function element_is_inside(Y::Matrix{Bool})
    m = map(s -> s-1, size(Y))
    E = zeros(Int8, m...)
    for i_2 in Base.axes(E, 2)
        for i_1 in Base.axes(E, 1)
            E[i_1, i_2] = sum(Y[i_1:i_1+1, i_2:i_2+1])
        end
    end  
    return E
end      

function element_is_inside(Y::Array{Bool,3})
    m = map(s -> s-1, size(Y))
    E = zeros(Int8, m...)
    for i_3 in Base.axes(E, 3)
        for i_2 in Base.axes(E, 2)
            for i_1 in Base.axes(E, 1)
                E[i_1, i_2, i_3] = sum(Y[i_1:i_1+1, i_2:i_2+1, i_3:i_3+1])
            end
        end  
    end
    return E
end      

"""
    active_splines(U::TensorProduct{2, <:SplineSpace}, E::Array)

Determine the active B-splines and boundary B-splines given a 
distance function ϕ.

A[k] == 0 => no physical elements in support of the function
A[k] == 1 => boundary function with a physical cut element in its support
A[k] == 2 => at least one physical element is fully inside support
A[k] == 3 => all physical elements are fully inside support
"""
function active_splines(S::SplineSpace, E::Vector{Int8})
    
    n = dimsplinespace(S)
    A = zeros(Int8, n...)
    
    supp = Support(S)

    for i in eachindex(A)
        s = supp[i]
        
        if any(E[s] .> 0)
            
            if any(E[s] .== 2) # inside
                
                A[i] = 2
                
                if all(E[s] .== 2)
                    A[i] = 3
                end
            
            else # boundary functions
                A[i] = 1
            end
        end
    end

    return A
end

active_splines(S::TensorProduct{1,<:SplineSpace}, E) = active_splines(S[1], E)

"""
    active_splines(U::TensorProduct{2, <:SplineSpace}, E::Array)

Determine the active B-splines and boundary B-splines given a 
distance function ϕ.

A[k] == 0 => no physical elements in support of the function
A[k] == 1 => boundary function with a physical cut element in its support
A[k] == 2 => at least one physical element is fully inside support
A[k] == 3 => all physical elements are fully inside support
"""
function active_splines(U::TensorProduct{2, <:SplineSpace}, E::Matrix{Int8})
    
    n = map(s -> dimsplinespace(s), U)
    A = zeros(Int8, n...)
    
    supp = map(s -> collect(Support(s)), U)

    for i_2 in Base.axes(A, 2)
        s2 = supp[2][i_2]
        
        for i_1 in Base.axes(A, 1)
            s1 = supp[1][i_1]
            
            if any(E[s1, s2] .> 0)
                
                if any(E[s1, s2] .== 4) # inside
                    
                    A[i_1, i_2] = 2
                    
                    if all(E[s1, s2] .== 4)
                        A[i_1, i_2] = 3
                    end
                
                else # boundary functions
                    A[i_1, i_2] = 1
                end
            end
        end
    end

    return A
end

"""
    active_splines(U::TensorProduct{Dim, <:SplineSpace}, E::Array)

Determine the active B-splines and boundary B-splines given a 
distance function ϕ.

A[k] == 0 => no physical elements in support of the function
A[k] == 1 => boundary function with a physical cut element in its support
A[k] == 2 => at least one physical element is fully inside support
A[k] == 3 => all physical elements are fully inside support
"""
function active_splines(U::TensorProduct{3, <:SplineSpace}, E::Array{Int8, 3})
    
    n = map(s -> dimsplinespace(s), U)
    A = zeros(Int8, n...)
    
    supp = map(s -> collect(Support(s)), U)

    for i_3 in Base.axes(A, 3)
        s3 = supp[3][i_3]

        for i_2 in Base.axes(A, 2)
            s2 = supp[2][i_2]
            
            for i_1 in Base.axes(A, 1)
                s1 = supp[1][i_1]
                
                if any(E[s1, s2, s3] .> 0)
                    
                    if any(E[s1, s2, s3] .== 8) # inside
                        
                        A[i_1, i_2, i_3] = 2
                        
                        if all(E[s1, s2, s3] .== 8)
                            A[i_1, i_2, i_3] = 3
                        end
                    
                    else # boundary functions
                        A[i_1, i_2, i_3] = 1
                    end
                end
            end
        end
    end
    return A
end

function active_splines(U::TensorProduct{Dim, <:SplineSpace}, phi) where Dim
    partition = CartesianProduct(u -> breakpoints(u), U)
    E = element_is_inside(phi, partition)
    return active_splines(U, E)
end

function active_splines(S::SplineSpace, phi)
    partition = breakpoints(S)
    E = element_is_inside(phi, partition)
    return active_splines(S, E)
end


"""
    find_extension_indices(active_functions, gperm, mi)

Find the extension indices that are `closest` to the boundary
function with CartesianIndex mi.
"""
function find_first_active_element(active_elements, gperm, mi)
    for ci in gperm(mi, k_ring)
        if active_elements[ci]
            return ci
        end
    end
end


"""
    find_extension_indices(active_functions, gperm, mi)

Find the extension indices that are `closest` to the boundary
function with CartesianIndex mi.
"""
function find_extension_indices(active_functions, gperm, mi)
    k_ring = 1
    while k_ring<20
        for ci in gperm(mi, k_ring)
            if has_black_neighbors_only(active_functions, ci)
                return ci
            end
        end
        k_ring+=1
    end
    throw(ArgumentError("Extension coefficients cannot be determined. Check if your input distance function is valid."))
end

function has_black_neighbors_only(active_functions, ci)
    for i in ci
        v = active_functions[i.I...]
        if v==1 || v==0
            return false
        end
    end
    return true
end