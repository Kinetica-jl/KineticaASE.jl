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


mutable struct FHIAimsBuilder
    calc_class::Py
    command::String
    xc::String
    species_dir::String
    maxiter::Int
    dispersion::String
    sc_accuracy_rho::Union{Nothing, Float64}
    sc_accuracy_forces::Union{Nothing, Float64}
    sc_accuracy_etot::Union{Nothing, Float64}
    sc_accuracy_eev::Union{Nothing, Float64}
end

function FHIAimsBuilder(;
        command::String="aims.x",
        xc::String="pbe",
        species_dir::String="./species_defaults/defaults_2020/tight",
        maxiter::Int=1000,
        dispersion::String="",
        sc_accuracy_rho::Union{Nothing, Float64}=nothing,
        sc_accuracy_forces::Union{Nothing, Float64}=nothing,
        sc_accuracy_etot::Union{Nothing, Float64}=nothing,
        sc_accuracy_eev::Union{Nothing, Float64}=nothing
    )
    if !isdir(species_dir)
        throw(ArgumentError("No species_dir found at $(species_dir)"))
    end

    aims = pyimport("ase.calculators.aims")
    return FHIAimsBuilder(aims.Aims, command, xc, species_dir, maxiter, dispersion,
                          sc_accuracy_rho, sc_accuracy_forces, sc_accuracy_etot, sc_accuracy_eev)
end

function (builder::FHIAimsBuilder)(dir::String, mult::Int, chg::Int, kwargs...)
    arg_dict = Dict(
        "aims_command" => builder.command,
        "outfilename" => joinpath(dir, "aims.out"),
        "xc" => builder.xc,
        "species_dir" => builder.species_dir,
        "sc_iter_limit" => string(builder.maxiter)
    )
    if !(builder.dispersion == "")
        if count(c->c==' ', builder.dispersion) == 0
            arg_dict[builder.dispersion] = ""
        else
            disptype, dispargs = split(builder.dispersion, limit=2)
            arg_dict[string(disptype)] = string(dispargs)
        end
    end
    if !isnothing(builder.sc_accuracy_forces)
        arg_dict["sc_accuracy_forces"] = builder.sc_accuracy_forces
    else
        arg_dict["compute_forces"] = ".true."
    end
    if !isnothing(builder.sc_accuracy_rho) arg_dict["sc_accuracy_rho"] = builder.sc_accuracy_rho end
    if !isnothing(builder.sc_accuracy_etot) arg_dict["sc_accuracy_etot"] = builder.sc_accuracy_etot end
    if !isnothing(builder.sc_accuracy_eev) arg_dict["sc_accuracy_eev"] = builder.sc_accuracy_eev end

    arg_dict["spin"] = mult > 1 ? "collinear" : "none"
    arg_dict["charge"] = chg
    
    return builder.calc_class(arg_dict...)
end
