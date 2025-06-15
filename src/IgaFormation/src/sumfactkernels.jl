"""
    sumfact!(A::Tuple{<:Array,<:Array,<:Array}, B, C, testfuns, trialfuns)

Apply sum factorization. `A` and `B` are allocated using `pre_allocate_arrays(...)`.
`testfuns` and `trialfuns` are tuples of arrays that represent the test and trial functions
evaluated at the quadrature nodes. It is assumed that the quadrature weights have been 
incorporated into the test functions
"""
function sumfact!(A::Tuple{<:DenseArray,<:DenseArray,<:DenseArray}, B, C, testfuns, trialfuns)
    sumfact_3d_stage3!(A[3], B[3], C, testfuns[3], trialfuns[3])
    sumfact_3d_stage2!(A[2], B[2], A[3], testfuns[2], trialfuns[2])
    sumfact_3d_stage1!(A[1], B[1], A[2], testfuns[1], trialfuns[1])
    return A[1]
end

function sumfact_3d_stage3!(A::AbstractArray{T,4}, B::AbstractMatrix{T}, C::AbstractArray{T,3}, testfuns::AbstractMatrix{T}, trialfuns::AbstractMatrix{T}) where T

    # check size input
    q_1, q_2, q_3 = size(C)
    @assert size(trialfuns,2) == size(testfuns,2) == q_3
    @assert size(B,1) == size(testfuns,1) && size(B,2) == size(trialfuns,1)
    @assert size(A) == (size(testfuns,1), size(trialfuns,1), q_1, q_2)

    # loop over quadrature points in 3-dir
    for k in Base.OneTo(q_3)

        # compute vector outer product
        form_outer_product!(B, testfuns, trialfuns, k)

        # add contribution point k to A
        ii = 1
        for j in Base.OneTo(q_2)
            for i in Base.OneTo(q_1)

                @inbounds c = C[i,j,k]
                ii = add_contribution_quadrature_point!(A, B, c, ii)

            end
        end

    end
end

function sumfact_3d_stage2!(A::AbstractArray{T,3}, B::AbstractMatrix{T}, C::AbstractArray{T,4}, testfuns::AbstractMatrix{T}, trialfuns::AbstractMatrix{T}) where T

    # check size input
    m, n, q_1, q_2 = size(C)
    @assert size(trialfuns,2) == size(testfuns,2) == q_2
    @assert size(B,1) == size(testfuns,1) && size(B,2) == size(trialfuns,1)
    @assert size(A) == (m*size(testfuns,1), n*size(trialfuns,1), q_1)

    # loop over quadrature points in 2-dir
    for j in Base.OneTo(q_2)

        # compute vector outer product at quadrature point j
        form_outer_product!(B, testfuns, trialfuns, j)

        # loop over indices of A
        ii = 1
        for i in Base.OneTo(q_1)

            for t in Base.OneTo(n)
                for q in axes(B,2)
                    for s in Base.OneTo(m)

                        @inbounds c = C[s,t,i,j]
                        ii = add_kron_contribution_quadrature_point!(A, B, c, ii, q)

                    end
                end
            end
        end

    end
end

function sumfact_3d_stage1!(A::AbstractMatrix{T}, B::AbstractMatrix{T}, C::AbstractArray{T,3}, testfuns::AbstractMatrix{T}, trialfuns::AbstractMatrix{T}) where T

    # check size input
    m, n, q_1 = size(C)
    @assert size(trialfuns,2) == size(testfuns,2) == q_1
    @assert size(B,1) == size(testfuns,1) && size(B,2) == size(trialfuns,1)
    @assert size(A) == (m*size(testfuns,1), n*size(trialfuns,1))

    # loop over quadrature points in 1-dir
    for i in Base.OneTo(q_1)

        # compute vector outer product at quadrture point i
        form_outer_product!(B, testfuns, trialfuns, i)

        # loop over indices of A
        ii = 1
        for t in Base.OneTo(n)
            for q in axes(B,2)
                for s in Base.OneTo(m)

                    @inbounds c = C[s,t,i]
                    ii = add_kron_contribution_quadrature_point!(A, B, c, ii, q)

                end
            end
        end
    end
end


"""
    sumfact!(A::Tuple{<:Array,<:Array,<:Array}, C, testfuns)

Apply sum factorization. `A` and `B` are allocated using `pre_allocate_arrays(...)`.
`testfuns` and `trialfuns` are tuples of arrays that represent the test and trial functions
evaluated at the quadrature nodes. It is assumed that the quadrature weights have been 
incorporated into the test functions.
"""
function sumfact!(A::Tuple{<:DenseArray,<:DenseArray,<:DenseArray}, C, testfuns)
    sumfact_3d_stage3!(A[3], C, testfuns[3])
    sumfact_3d_stage2!(A[2], A[3], testfuns[2])
    sumfact_3d_stage1!(A[1], A[2], testfuns[1])
    return A[1]
end

function sumfact_3d_stage3!(A::AbstractArray{T,3}, C::AbstractArray{T,3}, testfuns::AbstractMatrix{T}) where T

    # check size input
    q_1, q_2, q_3 = size(C)
    @assert size(testfuns,2) == q_3
    @assert size(A) == (size(testfuns,1), q_1, q_2)

    # loop over quadrature points in 3-dir
    for k in Base.OneTo(q_3)

        @inbounds B = view(testfuns, :, k)

        # add contribution point k to A
        ii = 1
        for j in Base.OneTo(q_2)
            for i in Base.OneTo(q_1)

                @inbounds c = C[i,j,k]
                ii = add_contribution_quadrature_point!(A, B, c, ii)

            end
        end

    end
end

function sumfact_3d_stage2!(A::AbstractArray{T,2}, C::AbstractArray{T,3}, testfuns::AbstractMatrix{T}) where T

    # check size input
    m, q_1, q_2 = size(C)
    @assert size(testfuns,2) == q_2
    @assert size(A) == (m*size(testfuns,1), q_1)

    # loop over quadrature points in 2-dir
    for j in Base.OneTo(q_2)

        @inbounds B = view(testfuns, :, j)

        # loop over indices of A
        ii = 1
        for i in Base.OneTo(q_1)

            for s in Base.OneTo(m)

                @inbounds c = C[s,i,j]
                ii = add_contribution_quadrature_point!(A, B, c, ii)

            end
        end

    end
end

function sumfact_3d_stage1!(A::AbstractVector{T}, C::AbstractArray{T,2}, testfuns::AbstractMatrix{T}) where T

    # check size input
    m, q_1 = size(C)
    @assert size(testfuns,2) == q_1
    @assert size(A) == (m*size(testfuns,1),)

    # loop over quadrature points in 1-dir
    for i in Base.OneTo(q_1)

        @inbounds B = view(testfuns, :, i)

        # loop over indices of A
        ii = 1
        for s in Base.OneTo(m)

            @inbounds c = C[s,i]
            ii = add_contribution_quadrature_point!(A, B, c, ii)

        end
    end
end

"""
    sumfact!(A::Tuple{<:Array,<:Array}, B, C, testfuns, trialfuns)

Apply sum factorization. `A` and `B` are allocated using `pre_allocate_arrays(...)`.
`testfuns` and `trialfuns` are tuples of arays that represent the test and trial functions
evaluated at the quadrature nodes. It is assumed that the quadrature weights have been 
incorporated into the test functions.
"""
function sumfact!(A::Tuple{<:DenseArray,<:DenseArray}, B, C, testfuns, trialfuns)
    sumfact_2d_stage2!(A[2], B[2], C, testfuns[2], trialfuns[2])
    sumfact_2d_stage1!(A[1], B[1], A[2], testfuns[1], trialfuns[1])
    return A[1]
end

function sumfact_2d_stage2!(A::AbstractArray{T,3}, B::AbstractMatrix{T}, C::AbstractArray{T,2}, testfuns::AbstractMatrix{T}, trialfuns::AbstractMatrix{T}) where T

    # check size input
    q_1, q_2 = size(C)
    @assert size(testfuns,2) == size(trialfuns,2) == q_2
    @assert size(B,1) == size(testfuns,1) && size(B,2) == size(trialfuns,1)
    @assert size(A) == (size(testfuns,1), size(trialfuns,1), q_1)

    # loop over quadrature points in 3-dir
    for j in Base.OneTo(q_2)

        # compute vector outer product
        form_outer_product!(B, testfuns, trialfuns, j)

        # add contribution point k to A
        ii = 1
        for i in Base.OneTo(q_1)

            @inbounds c = C[i,j]
            ii = add_contribution_quadrature_point!(A, B, c, ii)

        end

    end
end

function sumfact_2d_stage1!(A::AbstractMatrix{T}, B::AbstractMatrix{T}, C::AbstractArray{T,3}, testfuns::AbstractMatrix{T}, trialfuns::AbstractMatrix{T}) where T

    # check size input
    m, n, q_1 = size(C)
    @assert size(testfuns,2) == size(trialfuns,2) == q_1
    @assert size(B,1) == size(testfuns,1) && size(B,2) == size(trialfuns,1)
    @assert size(A) == (m*size(testfuns,1), n*size(trialfuns,1))

    # loop over quadrature points in 1-dir
    for i in Base.OneTo(q_1)

        # compute vector outer product at quadrture point i
        form_outer_product!(B, testfuns, trialfuns, i)

        # loop over indices of A
        ii = 1
        for t in Base.OneTo(n)
            for q in axes(B,2)
                for s in Base.OneTo(m)

                    @inbounds c = C[s,t,i]
                    ii = add_kron_contribution_quadrature_point!(A, B, c, ii, q)

                end
            end
        end
    end
end


"""
    sumfact!(A::Tuple{<:Array,<:Array}, C, testfuns)

Apply sum factorization. `A` is pre-allocated using `pre_allocate_rhs_arrays(...)`.
`testfuns` is a tuple of matrices that represent the test functions evaluated at 
the quadrature nodes. It is assumed that the quadrature weights have been incorporated
into the test functions.
"""
function sumfact!(A::Tuple{<:DenseArray,<:DenseArray}, C, testfuns)
    sumfact_2d_stage2!(A[2], C, testfuns[2])
    sumfact_2d_stage1!(A[1], A[2], testfuns[1])
    return A[1]
end

function sumfact_2d_stage2!(A::AbstractArray{T,2}, C::AbstractArray{T,2}, testfuns::AbstractMatrix{T}) where T

    # check size input
    q_1, q_2 = size(C)
    @assert size(testfuns,2) == q_2
    @assert size(A) == (size(testfuns,1), q_1)

    # loop over quadrature points in 3-dir
    for j in Base.OneTo(q_2)

        B = view(testfuns, :, j)

        # add contribution point k to A
        ii = 1
        for i in Base.OneTo(q_1)

            @inbounds c = C[i,j]
            ii = add_contribution_quadrature_point!(A, B, c, ii)

        end

    end
end

function sumfact_2d_stage1!(A::AbstractVector{T}, C::AbstractArray{T,2}, testfuns::AbstractMatrix{T}) where T

    # check size input
    m, q_1 = size(C)
    @assert size(testfuns,2)  == q_1
    @assert size(A,1) == m*size(testfuns,1)

    # loop over quadrature points in 1-dir
    for i in Base.OneTo(q_1)

        B = view(testfuns, :, i)

        # loop over indices of A
        ii = 1
        for s in Base.OneTo(m)

            @inbounds c = C[s,i]
            ii = add_contribution_quadrature_point!(A, B, c, ii)

        end

    end
end

# efficient evaluation of outer product of two vectors
@inline function form_outer_product!(B, a, b, k::Int)
    m=1
    for j in axes(B,2)
        @inbounds bj = b[j,k]
        for i in axes(B,1)
            @inbounds B[m] = a[i,k] * bj
            m+=1
        end
    end
end

# efficient evaluation of inner loop for one quadrature point
@inline function add_contribution_quadrature_point!(A, B, c, ii::Int)
    for q in axes(B,2)
        for p in axes(B,1)
            @inbounds A[ii] += c * B[p,q]
            ii+=1
        end
    end
    return ii
end

# efficient evaluation of inner loop with kron for one quadrature point
@inline function add_kron_contribution_quadrature_point!(A, B, c, ii::Int, q::Int)
    for p in axes(B,1)
        @inbounds A[ii] += c * B[p,q]
        ii+=1
    end
    return ii
end