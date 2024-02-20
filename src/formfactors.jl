# The part will be moved out of this module and should not be touched
function _pre_formfactor_grid(s::Integer)
    N = 3 * (s - 1)^2 + 9 * (s - 1) + 7

    sites = [SVector(0, 0)]

    edge1 = SVector(1, 0)
    edge2 = SVector(0, 1)
    edge3 = SVector(-1, 1)
    T1 = SVector(-1, 1)
    T2 = SVector(-1, 0)
    T3 = SVector(0, -1)

    for i in 1:s
        for j in 0:i-1
            vec = i * edge1 + T1 * j
            push!(sites, vec)
            push!(sites, -vec)
        end
        for j in 0:i-1
            vec = i * edge2 + T2 * j
            push!(sites, vec)
            push!(sites, -vec)
        end
        for j in 0:i-1
            vec = i * edge3 + T3 * j
            push!(sites, vec)
            push!(sites, -vec)
        end
    end

    sites
end

# This generates all form factors up to level (shell / 2) safe
function formfactor_grid(shells::Integer)
    _pre_formfactor_grid(2 * shells)
end

const ff_R1 = SVector(3.0/2.0,sqrt(3.0)/2.0)
const ff_R2 = SVector(3.0/2.0,-sqrt(3.0)/2.0)

########################
##### User part here ###
########################
"""
    filtered_formfactors(s)
Perform symmetry reduction on formfactors with `s` shells.
[`restore_formfactors!`](@restore_formfactors!) must undo this reduction!
"""
function filtered_formfactors(s)
    all_formfactors = formfactor_grid(s)

    filtered = [x for (k,x) in enumerate(all_formfactors) if k==1 || iseven(k)]
end

"""
    restore_formfactors!(all, reduced)
Reconstructs all elements in am array (with linear indexing) from the set where the precalculation has been carried out, i.e. symmetry redudancies are restored here.

# Arguments
- `all`: Stores the result. must be accessible via a single index
- `reduced`: Stores the input. Has a function type signature for the same index as `all`

# Examples
```julia-repl
julia> restore_formfactors!(sum_ph, i->precalcs[1, i, :] ./ Λ)
```
"""
function restore_formfactors!(all, reduced)
    for i in eachindex(all)
        if i==1
            all[i] = reduced(i)
        else
            map_i = (i÷2) + 1
            if i%2 == 0
                all[i] = reduced(map_i)
            else
                all[i] = conj(all[i-1])
            end
        end
    end

    nothing
end