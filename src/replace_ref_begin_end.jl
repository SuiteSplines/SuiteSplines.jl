"""
    replace_ref_begin_end!(ex)
Recursively replace occurrences of the symbols `:begin` and `:end` in a "ref" expression
(i.e. `A[...]`) `ex` with the appropriate function calls (`firstindex` or `lastindex`).
Replacement uses the closest enclosing ref, so
    A[B[end]]
should transform to
    A[B[lastindex(B)]]
"""
custom_replace_ref_begin_end!(ex) = custom_replace_ref_begin_end_!(ex, nothing)[1]
# replace_ref_begin_end_!(ex,withex) returns (new ex, whether withex was used)
function custom_replace_ref_begin_end_!(ex, withex)
    used_withex = false
    if isa(ex,Symbol)
        if ex === :begin
            withex === nothing && error("Invalid use of begin")
            return withex[1], true
        elseif ex === :end
            withex === nothing && error("Invalid use of end")
            return withex[2], true
        end
    elseif isa(ex,Expr)
        if ex.head === :ref
            ex.args[1], used_withex = custom_replace_ref_begin_end_!(ex.args[1], withex)
            S = isa(ex.args[1],Symbol) ? ex.args[1]::Symbol : gensym(:S) # temp var to cache ex.args[1] if needed
            used_S = false # whether we actually need S
            # new :ref, so redefine withex
            nargs = length(ex.args)-1
            if nargs == 0
                return ex, used_withex
            elseif nargs == 1
                # replace with lastindex(S)
                ex.args[2], used_S = custom_replace_ref_begin_end_!(ex.args[2], (:($firstindex($S)),:($lastindex($S))))
            else
                n = 1
                J = lastindex(ex.args)
                for j = 2:J
                    exj, used = custom_replace_ref_begin_end_!(ex.args[j], (:($firstindex($S,$n)),:($lastindex($S,$n))))
                    used_S |= used
                    ex.args[j] = exj
                    if isa(exj,Expr) && exj.head === :...
                        # splatted object
                        exjs = exj.args[1]
                        n = :($n + length($exjs))
                    elseif isa(n, Expr)
                        # previous expression splatted
                        n = :($n + 1)
                    else
                        # an integer
                        n += 1
                    end
                end
            end
            if used_S && S !== ex.args[1]
                S0 = ex.args[1]
                ex.args[1] = S
                ex = Expr(:let, :($S = $S0), ex)
            end
        else
            # recursive search
            for i = eachindex(ex.args)
                ex.args[i], used = custom_replace_ref_begin_end_!(ex.args[i], withex)
                used_withex |= used
            end
        end
    end
    ex, used_withex
end