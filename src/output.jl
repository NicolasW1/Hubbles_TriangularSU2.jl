const format_param_precision = 3
const format_momentum_precision = 6

function parameterString(p::t₁Params)
    st1 = format(p.t₁, precision=format_param_precision)
    smu = format(p.μ, precision=format_param_precision)

    "t1_" * st1 * "_" * "mu_" * smu
end
function parameterString(params::t₂Params)
    st2 = format(params.t₂, precision=format_param_precision)
    smu = format(params.μ, precision=format_param_precision)

    "t2_" * st2 * "_" * "mu_" * smu
end

function momentumString(q::SVector{2,T}) where {T}
    sqx = format(q[1], precision=format_momentum_precision)
    sqy = format(q[2], precision=format_momentum_precision)

    "qx_" * sqx * "_" * "qy_" * sqy
end

function outputFolder(extParams::ExternalInput)
    sparam = parameterString(extParams.params)

    if @isdefined(outputBasePath)
        path = joinpath(outputBasePath, sparam)
    else
        path = joinpath(@__DIR__, "tmp_output", sparam)
        @warn("Variable `outputBasePath` not defined! Output will be written to $(path). This is not safe and strongly disencouraged!")
    end

    mkpath(path)

    path
end

function outputFileName(bubble::Type, extParams::ExternalInput)
    smom = momentumString(extParams.qₑₓₜ)
    sBubble = string(bubble)

    sBubble * "_" * smom
end