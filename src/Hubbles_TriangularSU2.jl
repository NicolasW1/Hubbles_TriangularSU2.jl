module Hubbles_TriangularSU2

using StaticArrays, Formatting

include("will_be_another_module.jl")
export AbstractParameter, Bubble, TnFp, ∂¹aTnFp, ∂³aTnFp, ExternalInput

include("channels.jl")
include("parameter.jl")
include("formfactors.jl")

include("dispertion_relations.jl")
include("integrands.jl")

include("output.jl")

export ParticleParticle, ParticleHole
export t₁Params, t₂Params, Parameter, saveDict
export filtered_formfactors, restore_formfactors!
export bubble!, onsite_bubble

export parameterString, momentumString, outputFolder, outputFileName

end
