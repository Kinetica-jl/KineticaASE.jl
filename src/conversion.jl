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
    atoms_to_frame(atoms)

Converts an ASE Atoms object to an ExtXYZ frame.
"""
function atoms_to_frame(atoms::Py)
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
    return frame
end