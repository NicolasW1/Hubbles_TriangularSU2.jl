

# kernel = bubble without formfactor contribution
#####################################
############## Kernels ##############
#####################################

# the Δξ³ term vanishes by design, the threshold is set to (eps(Float64)^(1/3)/3!)
# consequently the next order contribution Δξ⁴/4! vanishes (exactly) for ∂ⁿaTnFp being order unity
const float_threshold = 0.000011003490984401003

# This threshold is the reason the two kernel functions are restricted to Float64, this number obviously changes depending on the Floating Point type

function kernel_Kplus(T::Float64, ξ₁::Float64, ξ₂::Float64)
    ξb = (ξ₁ + ξ₂) / (2 * T)
    Δξ = (ξ₁ - ξ₂) / (2 * T)

    if iszero(Δξ)
        κ = exp(-ξb)
        ∂¹aTnFp(ξb, κ) / T
    elseif abs(Δξ) < float_threshold
        κ = exp(-ξb)
        (∂¹aTnFp(ξb, κ) + Δξ^2/6 * ∂³aTnFp(ξb, κ)) / T
    else
        κ₁ = exp(-ξ₁/T)
        κ₂ = exp(-ξ₂/T)

        (TnFp(ξ₁/T, κ₁) - TnFp(ξ₂/T, κ₂)) / (ξ₁ - ξ₂)
    end
end

function kernel_Kminus(T::Float64, ξ₁::Float64, ξ₂::Float64)
    ξb = (ξ₁ + ξ₂) / (2 * T)
    Δξ = (ξ₁ - ξ₂) / (2 * T)

    if iszero(ξb)
        κ = exp(-Δξ)
        -∂¹aTnFp(Δξ, κ) / T
    elseif abs(ξb) < float_threshold
        κ = exp(-Δξ)
        -(∂¹aTnFp(Δξ, κ) + ξb^2/6 * ∂³aTnFp(Δξ, κ)) / T
    else
        κ₁ = exp(-ξ₁/T)
        κ₂ = exp(-ξ₂/T)

        -(TnFp(ξ₁/T, κ₁) + TnFp(ξ₂/T, κ₂)) / (ξ₁ + ξ₂)
    end
end

#####################################
############### evals ###############
#####################################

function evalKernel_PH(qₑₓₜ::SVector{2, T}, pₗₒₒₚ::SVector{2,T}, Λ::T, params::AbstractParameter{T}) where {T <: AbstractFloat}
    kernel_Kplus(Λ, ξ(pₗₒₒₚ + qₑₓₜ, params), ξ(pₗₒₒₚ, params))
end

function evalKernel_PP(qₑₓₜ::SVector{2, T}, pₗₒₒₚ::SVector{2,T}, Λ::T, params::AbstractParameter{T}) where {T <: AbstractFloat}
    kernel_Kminus(Λ, ξ(pₗₒₒₚ + qₑₓₜ, params), ξ(-pₗₒₒₚ, params))
end

#####################################
############## Bubbles ##############
#####################################

function bubble!(result, ::Type{ParticleParticle}, pₗₒₒₚ::SVector{2,T}, ext::ExternalInput{T, I, P}) where {T, I, P}
    kernel = evalKernel_PP(ext.qₑₓₜ, pₗₒₒₚ, ext.Λ, ext.params)
    
    for i ∈ eachindex(result)
        result[i] = kernel * formfactor(pₗₒₒₚ, ext.ffs[i])
    end

    nothing
end

function onsite_bubble(::Type{ParticleParticle}, pₗₒₒₚ::SVector{2,T}, ext::ExternalInput{T, I, P}) where {T, I, P}
    evalKernel_PP(ext.qₑₓₜ, pₗₒₒₚ, ext.Λ, ext.params)
end

function bubble!(result, ::Type{ParticleHole}, pₗₒₒₚ::SVector{2,T}, ext::ExternalInput{T, I, P}) where {T, I, P}
    kernel = evalKernel_PH(ext.qₑₓₜ, pₗₒₒₚ, ext.Λ, ext.params)

    for i ∈ eachindex(result)
        result[i] = kernel * formfactor(pₗₒₒₚ, ext.ffs[i])
    end

    nothing
end
function onsite_bubble(::Type{ParticleHole}, pₗₒₒₚ::SVector{2,T}, ext::ExternalInput{T, I, P}) where {T, I, P}
    evalKernel_PH(ext.qₑₓₜ, pₗₒₒₚ, ext.Λ, ext.params)
end