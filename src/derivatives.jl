# overloading of AbstractMappings.evalkernel! which enables computation
# on CartesianProduct grids of points using the @evaluate! macro.
for Dim in 1:3
    @eval function AbstractMappings.evalkernel_grad_imp!(op::Val{:(=)}, y::EvaluationSet{$Dim,1}, x, nurbs::Nurbs{$Dim})
        
        # evaluate Nurbs function
        S, W = nurbs.splinefun, nurbs.weightfun
        @evaluate s  = S(x)
        @evaluate w  = W(x)
        _eval_nurbs!(s, w)

        # evaluate Nurbs gradient
        for dir in 1:dimension(nurbs)
            dS, dW = partial_derivative(S, 1, dir), partial_derivative(W, 1, dir)
            ds = y.data[dir]
            @evaluate! ds  = dS(x)
            @evaluate dw  = dW(x)
            _eval_nurbs_derivative!(s, ds, w, dw)
        end
    end

    # overloading of AbstractMappings.evalkernel_hessian_imp! which enables computation
    # on CartesianProduct grids of points using the @evaluate and @evaluate! macro.
    @eval function AbstractMappings.evalkernel_hessian_imp!(op::Val{:(=)}, Δs::EvaluationSet{$Dim,$Dim}, x, nurbs::Nurbs{$Dim})

        # evaluate Nurbs function
        S, W = nurbs.splinefun, nurbs.weightfun
        @evaluate s = S(x)
        @evaluate w = W(x)
        _eval_nurbs!(s, w)

        # evaluate gradients
        @evaluate ∇s = Gradient(S)(x)
        @evaluate ∇w = Gradient(W)(x)
        _eval_nurbs_gradient!(s, ∇s, w, ∇w)

        # evaluate hessian
        @evaluate! Δs = Hessian(S)(x)
        @evaluate Δw = Hessian(W)(x)
        _eval_nurbs_hessian!(s, ∇s, Δs, w, ∇w, Δw)
    end
end

# evaluate nurbs derivatives inplace
function _eval_nurbs_derivative!(Y, dY, W, dW)
    dY .-= dW .* Y
    dY ./= W
end

# evaluate nurbs gradient inplace
function _eval_nurbs_gradient!(Y, ∇Y, W, ∇W)
    for k in 1:length(∇Y.data)
        _eval_nurbs_derivative!(Y, ∇Y.data[k], W, ∇W.data[k])
    end
end


# evaluate nurbs hessian inplace
function _eval_nurbs_hessian!(Y, ∇Y, ΔY, W, ∇W, ΔW)
    for j in 1:size(ΔY.data,2)
        for i in 1:size(ΔY.data,1)
            ΔY.data[i,j] .-= ΔW.data[i,j] .* Y
            ΔY.data[i,j] .-= ∇W.data[i] .* ∇Y.data[j]
            ΔY.data[i,j] .-= ∇W.data[j] .* ∇Y.data[i]
            ΔY.data[i,j] ./= W
        end
    end
end
