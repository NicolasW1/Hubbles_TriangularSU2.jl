struct Parameter{T<:AbstractFloat} <: AbstractParameter{T}
    t₁::T
    t₂::T
    t₃::T
    μ::T
end
saveDict(p::Parameter) = Dict("t1" => p.t₁, "t2" => p.t₂, "t3" => p.t₃, "mu" => p.μ)

struct t₁Params{T<:AbstractFloat} <: AbstractParameter{T}
    t₁::T
    μ::T
end
Parameter(p::t₁Params{T}) where {T} = Parameter(p.t₁, zero(T), zero(T), p.μ)
t₁Params(μ::AbstractFloat) = t₁Params(one(μ), μ)
t₁Params(μ::Number) = t₁Params(float(μ))
saveDict(p::t₁Params) = saveDict(Parameter(p))

struct t₂Params{T<:AbstractFloat} <: AbstractParameter{T}
    t₁::T
    t₂::T
    μ::T
end
Parameter(p::t₂Params{T}) where {T} = Parameter(p.t₁, p.t₂, zero(T), p.μ)
t₂Params(t₂::AbstractFloat) = t₂Params(one(t₂), t₂, 2 * (t₂ + 1))
t₂Params(t₂::Number) = t₂Params(float(t₂))
saveDict(p::t₂Params) = saveDict(Parameter(p))

