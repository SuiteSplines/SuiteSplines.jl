export @evaluate, @evaluate!

const EvaluationCache = LRU{Tuple{AbstractMapping,Int},EvaluationSet}(maxsize=1000000, by=sizeof)

function get_evaluation_cache(mapping::AbstractMapping, x)
    maxlength = base_2_ceil(Int, prod(size(x)))
    key = (mapping, maxlength)
    
    # add key-value pair to cache if needed
    if !haskey(EvaluationCache, key)
        T, Dim = eltype(x[1]), ndims(x)
        S1, S2 = codimension(mapping)
        y = ntuple(k -> ResizableArray{T,Dim}(undef, capacity=maxlength), S1*S2)
        push!(EvaluationCache, key => EvaluationSet{S1,S2}(y...))
    end

    # extract cached array
    y = EvaluationCache[key]

    # resize result and return
    [resize!(y.data[k], size(x)) for k in 1:length(y.data)]
    return y
end

Base.size(x::NTuple{Dim,<:AbstractVector}) where Dim = ntuple(k->length(x[k]), Dim)

# kernels that needs to be implemented by every sub-type of AbstractMapping in order 
# to use the @evaluation and @evaluation! macro's
function evalkernel! end

# kernels that needs to be implemented by every sub-type of Gradient in order 
# to use the @evaluation and @evaluation! macro's
function evalkernel_grad_imp! end

# kernels that needs to be implemented by every sub-type of Hessian in order 
# to use the @evaluation and @evaluation! macro's
function evalkernel_hessian_imp! end

"""
    @evaluate Y  = f(x)

Fast evaluation routine for evaluating a mapping on a grid  ``x ...``. The
result is stored in array ``Y``. Allocation is done automatically and EvaluationSet's
of the correct size and type are cached for each mapping using an LRUCache.
"""
macro evaluate(ex)
    local op = esc(Val(ex.head))
    local y = esc(ex.args[1])
    local f = esc(ex.args[2].args[1])
    local x = esc(ex.args[2].args[2])
    @assert ex.head==:(=) "The signature should be 'y = f(x)'."
    quote
        $y = evaluate_imp($x, $f)
    end
end

"""
    @evaluate! Y = f(x)
    @evaluate! Y += f(x)
    @evaluate! Y -= f(x)
    @evaluate! Y *= f(x)
    @evaluate! Y /= f(x)

Fast update routine for evaluating a mapping on a grid  ``x ...``. array
`Y` needs to be preallocated.
"""
macro evaluate!(ex)
    local op = esc(Val(ex.head))
    local y = esc(ex.args[1])
    local f = esc(ex.args[2].args[1])
    local x = esc(ex.args[2].args[2])
    @assert ex.head==:(=) || ex.head==:(+=) || ex.head==:(-=) || ex.head==:(*=) || ex.head==:(/=) "The signature should be 'y = f(x)', 'y += f(x)', 'y -= f(x)', 'y *= f(x)' or 'y /= f(x)'."
    quote
        evaluate_imp!($op, $y, $x, $f)
    end
end

@inline function evaluate_imp(x, f)
    y = get_evaluation_cache(f, x)
    evalkernel!(Val(:(=)), y, x, f)
    return y
end

@inline function evaluate_imp!(op, y::EvaluationSet, x, f)
    evalkernel!(op, y, x, f)
    return y
end

@inline function evaluate_imp!(op, y, x, f)
    S1, S2 = codimension(f)
    evalkernel!(op, EvaluationSet{S1,S2}(y), x, f)
    return y
end

@inline function evaluate_imp(x, f::Pairing)
    @evaluate y = f.mapping(x)
    @evaluate z = f.field(y)
    return z
end

@inline function evaluate_imp!(op, z::EvaluationSet, x, f::Pairing)
    @evaluate y = f.mapping(x)
    evaluate_imp!(op, z, y, f.field)
    return z
end

# Specialized evaluation routines for Normal to a 3D surface
function evalkernel!(::Val{:(=)}, y::EvaluationSet{3,1}, x, n::Normal{2,3})
    @evaluate ∇f = Gradient(n.mapping)(x)
    for k in 1:length(x)
        y[k] = cross(∇f[k][1,:], ∇f[k][2,:])
    end
    return y
end

# Specialized evaluation routines for Normal to a 2D curve
function evalkernel!(::Val{:(=)}, y::EvaluationSet{2,1}, x, n::Normal{1,2})
    @evaluate ∇f = Gradient(n.mapping)(x)
    for k in 1:length(x)
        y[k] = [∇f[k][2]; -∇f[k][1]]
    end
    return y
end

# Specialized evaluation routines for Normal to a 3D surface
function evalkernel!(::Val{:(=)}, y::EvaluationSet{3,1}, x, n::UnitNormal{2,3})
    @evaluate! y = Normal(n.mapping)(x)
    for k in 1:length(x)
        y[k] /= norm(y[k]) 
    end
    return y
end

# Specialized evaluation routines for Normal to a 2D curve
function evalkernel!(::Val{:(=)}, y::EvaluationSet{2,1}, x, n::UnitNormal{1,2})
    @evaluate! y = Normal(n.mapping)(x)
    for k in 1:length(x)
        y[k] /= norm(y[k]) 
    end
    return y
end

# Specialized evaluation routines for Vol - case Dim == Codim
for Dim in 1:3
    @eval function evalkernel!(::Val{:(=)}, y::EvaluationSet{1,1}, x, J::Vol{$Dim,<:GeometricMapping{$Dim,$Dim}})
        @evaluate ∇f = Gradient(J.mapping)(x)
        for k in 1:length(x)
            y[k] = det(∇f[k]) 
        end
        return y
    end
end

# Specialized evaluation routines for Vol - case Dim != Codim
function evalkernel!(::Val{:(=)}, y::EvaluationSet{1,1}, x, J::Vol{1,<:GeometricMapping{1,2}})
    @evaluate n = Normal(J.mapping)(x)
    for k in 1:length(x)
        y[k] = norm(n[k])
    end
    return y
end
function evalkernel!(::Val{:(=)}, y::EvaluationSet{1,1}, x, J::Vol{2,<:GeometricMapping{2,3}})
    @evaluate n = Normal(J.mapping)(x)
    for k in 1:length(x)
        y[k] = norm(n[k])
    end
    return y
end


using Base.Cartesian
for op in [:(=), :(+=), :(-=), :(*=), :(/=)]
    
    local S = Val{op}

    # Below are a set of general routines that are called via the macro @evaluate and @evaluate!
    # Overload evalkernel!, evalkernel_grad_imp!, evalkernel_hessian_imp! for your particular
    # datastructures in order to use the evaluation macro's. This is already done for standard
    # functions below

    # fast evaluation of an AbstractMapping on an arbitrary set of points
    for Dim in 1:3

        @eval function evalkernel!(op::$S, y::EvaluationSet{1,1}, x, mapping::ScalarMapping{$Dim})
            evalkernel_imp!(op, y.data[1], x, mapping)
        end

        for S2 in 1:3
            for S1 in 1:3
                @eval function evalkernel!(op::$S, y::EvaluationSet{$S1,$S2}, x, mapping::AbstractMapping{$Dim,$S1,$S2})
                    for k in 1:$S1*$S2
                        evalkernel_imp!(op, y.data[k], x, mapping[k])
                    end
                end
            end
        end
    end

    # fast evaluation of Gradients of mappings
    for Dim in 1:3

        @eval function evalkernel!(op::$S, Y::EvaluationSet{$Dim,1}, X, f::Gradient{$Dim,1,<:ScalarMapping{$Dim}})
            evalkernel_grad_imp!(op, Y, X, f.mapping)
            return Y
        end

        for Codim in 1:3
            @eval function evalkernel!(op::$S, Y::EvaluationSet{$Dim,$Codim}, X, f::Gradient{$Dim,$Codim})
                for k in 1:$Codim
                    y = EvaluationSet{$Dim,1}(Y.data[:,k]...)
                    evalkernel_grad_imp!(op, y, X, f.mapping[k])
                end
                return Y
            end
        end
    end

    # fast evaluation of the Gradient of a mapping
    for Dim in 1:3

        @eval function evalkernel!(op::$S, Y::EvaluationSet{$Dim,$Dim}, X, f::Hessian{$Dim,<:ScalarMapping{$Dim}})
            evalkernel_hessian_imp!(op, Y, X, f.mapping)
            return Y
        end

        @eval function evalkernel!(op::$S, Y::EvaluationSet{$Dim,$Dim}, X, f::Hessian{$Dim})
            evalkernel_hessian_imp!(op, Y, X, f.mapping.data[1])
            return Y
        end
    end

    # Fall-back methods for fast evaluation of the Gradient and Hessian of a mapping
    for Dim in 1:3

        # fall-back method that evaluates component-wise
        @eval function evalkernel_grad_imp!(op::$S, Y::EvaluationSet{$Dim,1}, X, f)
            for k in 1:$Dim
                evalkernel_imp!(op, Y.data[k], X, IgaBase.partial_derivative(f, 1, k))
            end
        end

        # fall-back method that evaluates component-wise
        @eval function evalkernel_hessian_imp!(op::$S, Y::EvaluationSet{$Dim,$Dim}, X, f)
            for j in 1:$Dim
                for i in 1:$Dim
                    ders = ntuple(l -> Int(l==i)+Int(l==j), $Dim)
                    evalkernel_imp!(op, Y.data[i,j], X, IgaBase.partial_derivative(f, ders))
                end
            end
        end
    end

    # evalkernel!, evalkernel_grad_imp!, evalkernel_hessian_imp! are overloaded for regular
    # functions in order to use the @evaluation and @evaluation! macro.

    # Function evaluation in 1D
    local kernel_expression = Expr(op, Expr(:ref, :Y, :k), Expr(:call, :func, :x))
    @eval function evalkernel_imp!(::$S, Y::AbstractVector, X::AbstractVector, func::ScalarFunction)
        @assert length(Y)==length(X)
        for k in 1:length(X)
            x = X[k]
            @inbounds $kernel_expression
        end
        return Y
    end

    for Dim in 1:3

        # fast evaluation of standard Julia functions on a Cartesian grid of points
        # defined by a tuple or CartesianProduct of vectors
        local kernel_expression = Expr(op, :(@nref $Dim Y i), :(@ncall $Dim func x))
        @eval function evalkernel_imp!(::$S, Y::AbstractArray{T,$Dim}, X::CartesianProduct{$Dim}, func::ScalarFunction) where T
            @assert size(Y)==size(X)
            @nloops $Dim i Y d -> x_d = X.data[d][i_d] begin
                @inbounds $kernel_expression
            end
            return Y
        end

        # fast evaluation of standard Julia functions on an Arbitrary set of points
        # defined by an EvaluationSet
        for Codim in 1:3
            local kernel_expression = Expr(op, :(@nref $Codim Y i), :(@ncall $Dim f x))
            @eval function evalkernel_imp!(::$S, Y::AbstractArray{T,$Codim}, X::EvaluationSet{1,$Dim}, f::ScalarFunction) where T
                @assert size(Y)==size(X)
                @nloops $Codim i Y begin
                    y = @nref $Codim X i
                    x = @nextract $Dim x y
                    @inbounds $kernel_expression
                end
                return Y
            end
        end

        # fast evaluation of Gradients of functions on a CartesianProduct
        local kernel_expression = Expr(op, :(@nref $Dim A i), :(save[l]))
        @eval function evalkernel_grad_imp!(op::$S, Y::EvaluationSet{$Dim,1}, X, f::ScalarFunction{$Dim})
            g = (x) -> ForwardDiff.gradient((x) -> f(x...), x)
            @nloops $Dim i X begin
                x = @nref $Dim X i
                save = g(SVector(x))
                for l in 1:length(save)
                    A = Y.data[l]
                    $kernel_expression
                end
            end
        end

        # fast evaluation of Hessians of functions on a CartesianProduct
        local kernel_expression = Expr(op, :(@nref $Dim A i), :(save[l]))
        @eval function evalkernel_hessian_imp!(::$S, Y::EvaluationSet{$Dim,$Dim}, X::CartesianProduct{$Dim}, f::ScalarFunction)
            h = (x) -> ForwardDiff.hessian((x) -> f(x...), x)
            @nloops $Dim i X begin
                x = @nref $Dim X i
                save = h(SVector(x))
                for l in 1:length(save)
                    A = Y.data[l]
                    $kernel_expression
                end
            end
            return Y
        end

    end # for Dim in 1:3

end # for op in [:(=), :(+=), :(-=), :(*=), :(/=)]




# get the closest number higher number that is a base of 2. 
function base_2_ceil(T::Type{<:Real}, a)
    return 2^ceil(T,log(2,a))
end
