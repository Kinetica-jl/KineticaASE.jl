mutable struct EMTBuilder
    calc_class::Py
    fixed_cutoff::Bool
end

function EMTBuilder()
    emt = pyimport("ase.calculators.emt")
    return EMTBuilder(emt.EMT, true)
end

function (builder::EMTBuilder)(dir::String, mult::Int, chg::Int, kwargs...)
    return builder.calc_class(fixed_cutoff=builder.fixed_cutoff)
end



mutable struct NWChemDFTBuilder
    calc_class::Py
    command::String
    xc::String
    basis::Union{String, Dict{String, String}}
    maxiter::Int
    convergence::String
    adft::Bool
    memory::String
end

function NWChemDFTBuilder(; 
        command::String="nwchem PREFIX.nwi > PREFIX.nwo",
        xc::String="becke97", 
        basis::Union{String, Dict{String, String}}="3-21G",
        maxiter::Int=50,
        convergence::String="",
        adft::Bool=true, 
        memory::String="1024 mb"
    )
    nwchem = pyimport("ase.calculators.nwchem")
    return NWChemDFTBuilder(nwchem.NWChem, command, xc, basis, maxiter, convergence, adft, memory)
end

function (builder::NWChemDFTBuilder)(dir::String, mult::Int, chg::Int, kwargs...)
    dft_dict = Dict(
        "xc" => builder.xc,
        "mult" => mult,
        "maxiter" => builder.maxiter
    )
    if builder.adft dft_dict["adft"] = nothing end
    if !(builder.convergence == "") dft_dict["convergence"] = builder.convergence end

    calc = builder.calc_class(
        memory=builder.memory,
        dft=pydict(dft_dict),
        basis=builder.basis
    )
    calc.command = builder.command
    return calc
end