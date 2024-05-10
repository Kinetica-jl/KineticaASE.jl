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
    xc::String
    basis::Union{String, Dict{String, String}}
    adft::Bool
    memory::String
end

function NWChemDFTBuilder(; 
        xc::String="pbe", 
        basis::Union{String, Dict{String, String}}="3-21G",
        adft::Bool=true, 
        memory::String="1024mb"
    )
    nwchem = pyimport("ase.calculators.nwchem")
    return NWChemDFTBuilder(nwchem.NWChem, xc, basis, adft, memory)
end

function (builder::NWChemDFTBuilder)(dir::String, mult::Int, chg::Int, kwargs...)
    dft_dict = Dict(
        "xc" => builder.xc,
        "mult" => mult
    )
    if builder.adft dft_dict["adft"] = nothing end

    return builder.calc_class(
        memory=builder.memory,
        dft=dft_dict,
        basis=builder.basis
    )
end