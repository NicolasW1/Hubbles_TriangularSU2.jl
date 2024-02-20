function εξ₁(q::SVector{2,T}) where {T}
    2 * (2 * cos((3 * q[1]) / 2) * cos((sqrt(T(3)) * q[2]) / 2) + cos(sqrt(T(3)) * q[2]))
end
function εξ₂(q::SVector{2,T}) where {T}
    2 * (cos(3 * q[1]) + 2 * cos((3 * q[1]) / 2) * cos((3 * sqrt(T(3)) * q[2]) / 2))
end
function εξ₃(q::SVector{2,T}) where {T}
    2 * (2 * cos(3 * q[1]) * cos(sqrt(T(3)) * q[2]) + cos(2 * sqrt(T(3)) * q[2]))
end

function ξ(q::SVector{2,T}, params::Parameter{TP}) where {T,TP}
    -params.t₁ * εξ₁(q) - params.t₂ * εξ₂(q) - params.t₃ * εξ₃(q) - params.μ
end

# specalized versions for parameter sets where only a subset of possible interactions in non-zero
# this is for performance
function ξ(q::SVector{2,T}, params::t₁Params{TP}) where {T,TP}
    -params.t₁ * εξ₁(q) - params.μ
end
function ξ(q::SVector{2,T}, params::t₂Params{TP}) where {T,TP}
    -params.t₁ * εξ₁(q) - params.t₂ * εξ₂(q) - params.μ
end