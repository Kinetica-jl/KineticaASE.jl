"""
    frame_to_atoms(frame[, charges=nothing, magmoms=nothing])

Converts an ExtXYZ frame to an ASE Atoms object.

If the optional argument `charges` is provided with an 
integer array of formal charges on the atoms in `frame`,
these will be applied to the resulting Atoms object through
`Atoms.set_initial_charges()`.

If the optional argument `magmoms` is provided with an 
array of initial magnetic moments for the atoms in `frame`,
these will be applied to the resulting Atoms object through
`Atoms.set_initial_magnetic_moments()`.
"""
function frame_to_atoms(frame::Dict{String, Any}, charges=nothing, magmoms=nothing)
    symbols = join(frame["arrays"]["species"])
    positions = np.asarray(frame["arrays"]["pos"]')
    atoms = ase.Atoms(symbols, positions=positions)
    if !isnothing(charges)
        if pylen(atoms) != length(charges)
            throw(ArgumentError("Number of formal charges in `charges` must match number of atoms."))
        else
            atoms.set_initial_charges(charges)
        end
    end
    if !isnothing(magmoms)
        if pylen(atoms) != length(magmoms)
            throw(ArgumentError("Number of initial magnetic moments in `magmoms` must match number of atoms."))
        else
            atoms.set_initial_magnetic_moments(np.asarray(magmoms))
        end
    end
    atoms.info = frame["info"]
    return atoms
end


"""
    atoms_to_frame(atoms[, ase_energy])

Converts an ASE Atoms object to an ExtXYZ frame.

Optionally accepts values `ase_energy and `inertias`
that can fill the "energy_ASE" and "inertias" keys
of the resulting frame's "info" Dict. `ase_energy`
should be a Julia float in eV and `inertias` should
be a `Vector{Float64}` for proper compatibility.
"""
function atoms_to_frame(atoms::Py, ase_energy=nothing, inertias=nothing)
    symbols = pyconvert(Vector{String}, atoms.get_chemical_symbols())
    positions = pyconvert(Matrix, atoms.get_positions().T)
    frame = Dict{String, Any}(
        "N_atoms" => length(symbols),
        "arrays" => Dict{String, Any}(
            "species" => symbols,
            "pos" => positions
        ),
        "info" => pyconvert(Dict{String, Any}, atoms.info)
    )
    if !isnothing(ase_energy) frame["info"]["energy_ASE"] = ase_energy end
    if !isnothing(inertias) frame["arrays"]["inertias"] = inertias end
    return frame
end


"""
"""
function get_reverse_rhash(sd::SpeciesData, rd::RxData, rid)
    reacs = []
    for (i, sid) in enumerate(rd.id_reacs[rid])
        for _ in 1:rd.stoic_reacs[rid][i] 
            push!(reacs, sd.toStr[sid])
        end
    end
    sort!(reacs)
    prods = []
    for (i, sid) in enumerate(rd.id_prods[rid])
        for _ in 1:rd.stoic_prods[rid][i] 
            push!(prods, sd.toStr[sid])
        end
    end
    sort!(prods)

    forw_rhash = stable_hash(vcat(reacs, prods))
    @assert rd.rhash[rid] == forw_rhash
    rev_rhash = stable_hash(vcat(prods, reacs))
    return rev_rhash
end


"""
"""
function imaginary_ve_tol(imaginary_freq_tol)
    return (imaginary_freq_tol^-0.5) * ASEConstants.hbar * ASEConstants.m / sqrt(ASEConstants.e * ASEConstants.amu)
end