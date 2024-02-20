module Hubbles_TriangularSU2

using StaticArrays, Formatting

include("will_be_another_module.jl")
export AbstractParameter, Bubble, TnFp, ∂¹aTnFp, ∂³aTnFp, ExternalInput

include("parameter.jl")
include("channels.jl")
include("dispertion_relations.jl")

include("output.jl")

export ParticleParticle, ParticleHole
export t₁Params, t₂Params, Parameter, saveDict
export parameterString, momentumString, outputFolder, outputFileName

end
