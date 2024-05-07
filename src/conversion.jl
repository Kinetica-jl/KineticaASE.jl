"""
    frame_to_atoms(frame)

Converts an ExtXYZ frame to an ASE Atoms object.
"""
function frame_to_atoms(frame::Dict{String, Any})
    symbols = join(frame["arrays"]["species"])
    positions = frame["arrays"]["pos"]'
    atoms = ase.Atoms(symbols, positions=positions)
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
            "species" => species,
            "pos" => positions
        ),
        "info" => pyconvert(Dict{String, Any}, atoms.info)
    )
    if !isnothing(ase_energy) frame["info"]["energy_ASE"] = ase_energy end
    if !isnothing(inertias) frame["info"]["inertias"] = inertias end
    return frame
end