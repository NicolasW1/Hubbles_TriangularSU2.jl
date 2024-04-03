module Hubbles_TriangularSU2

using StaticArrays, Formatting
using HubblesPrerequisites

include("channels.jl")
include("parameter.jl")
include("formfactors.jl")

include("dispertion_relations.jl")
include("integrands.jl")

include("output.jl")

export ParticleParticle, ParticleHole
export t₁Params, t₂Params, Parameter
export filtered_formfactors, restore_formfactors!
export bubble!, onsite_bubble

export saveDict, parameterString, momentumString, outputFolder, outputFileName

struct Model{T1, T2, F1, F2, F3, F4, FM1, FM2, FM3, FM4, FM5} <: AbstractModel
    bubbles::T1
    Parameter::T2

    filtered_formfactors::F1
    restore_formfactors!::F2
    bubble!::F3
    onsite_bubble::F4

    saveDict::FM1
    parameterString::FM2
    momentumString::FM3
    outputFolder::FM4
    outputFileName::FM5
end

model = Model((ParticleParticle, ParticleHole), Parameter, filtered_formfactors, restore_formfactors!, bubble!, onsite_bubble, saveDict, parameterString, momentumString, outputFolder, outputFileName)

export Model, model

end
