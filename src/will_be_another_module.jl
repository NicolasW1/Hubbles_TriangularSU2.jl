# struct ExternalInput{T,X,P}
#     qₑₓₜ::SVector{2,T}
#     ffs::X
#     Λ::T
#     params::P
# end
# ExternalInput(ext::ExternalInput{T,X,P}, Λ::T) where {T,X,P} = ExternalInput(ext.qₑₓₜ, ext.ffs, Λ, ext.params)