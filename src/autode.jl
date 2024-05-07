
"""
    autode_to_frame(ademol::Py[, info_dict=nothing])

Converts an autodE `Molecule` to an ExtXYZ frame.

Reads atoms from `ademol` and constructs an ExtXYZ frame with
them. Frame defaults to having an empty `info` field, although
this can be populated by passing a `Dict{String, Any}` to the
`info_dict` keyword argument.
"""
function autode_to_frame(ademol::Py; info_dict=nothing)
    na = pylen(ademol.atoms)
    symbols = [pyconvert(String, ademol.atoms[i].atomic_symbol) for i in 0:na-1]
    coords = reduce(hcat, [
        [
            pyconvert(Float64, mol.atoms[i].coord.x),
            pyconvert(Float64, mol.atoms[i].coord.y),
            pyconvert(Float64, mol.atoms[i].coord.z)
        ] for i in 0:na-1
    ])
    
    info = isnothing(info_dict) ? Dict{String, Any}() : info_dict
    arrays = Dict{Sting, Any}("species" => symbols, "pos" => coords)
    frame = Dict{String, Any}("N_atoms" => na, "info" => info, "arrays" => arrays)
    return frame
end

"""
    frame_to_autode(frame::Dict{String, Any})

Converts an ExtXYZ frame to an autodE `Molecule`.

Writes a temporary xyz file from `frame` and reads it back in
with autodE. There doesn't seem to be a way around writing to
disk here, as autodE only detects an xyz input when reading a
file ending in '.xyz', which can't be replicated with an in-
memory IO buffer.
"""
function frame_to_autode(frame::Dict{String, Any})
    f = joinpath(tempdir(), "frame.xyz")
    write_frame(f, frame)
    mol = ade.Molecule(f)
    rm(f)
    return mol
end

"""
    autode_conformer_search!(sd::SpeciesData, i::Int)

Performs a conformer search for the species at `sd.xyz[i]`, updating its geometry.

Constructs an initial guess of species geometry from its SMILES,
then runs an autodE conformer search with xTB as the energetic
driver. Finds the lowest energy conformer and writes it back
to `sd.xyz[i]`.

Also populates `sd.cache` with autodE-derived values for symmetry
number, spin and geometry for later use in TST calculations.
"""
function autode_conformer_search!(sd, i)
    mol = ade.Molecule(smiles=sd.toStr[i])
    if sd.xyz[i]["N_atoms"] > 2
        @debug "Searching for conformers of species $i: $(sd.toStr[i])"
        mol.find_lowest_energy_conformer()
        n_confs_found = pylen(mol.conformers)
        @debug "$(n_confs_found) conformers found."

        @debug "Writing lowest energy conformer to SpeciesData."
        sd.xyz[i] = autode_to_frame(mol.conformers(lowest_e_id-1); info_dict=sd.xyz[i]["info"])
        sd.xyz[i]["info"]["energy"] = pyconvert(Float64, mol.energy.to("ev").real)
    else
        @debug "Species $i ($(sd.toStr[i])) too small for conformer search, skipping."
    end
    
    sd.cache[:symmetry][i] = pyconvert(Int, mol.symmetry_number)
    sd.cache[:spin][i] = (pyconvert(Int, mol.mult)-1)/2
    if pyconvert(Int, mol.n_atoms) == 1
        sd.cache[:geometry][i] = 0
    else
        if pyconvert(Bool, mol.is_linear())
            sd.cache[:geometry][i] = 1
        else
            sd.cache[:geometry][i] = 2
        end
    end

    rm(joinpath(pwd(), "conformers"), recursive=true)
end

"""
    autode_NCI_conformer_search(frames::Vector{Dict{String, Any}})
    
Performs a search for the lowest energy non-covalent interacting reaction complex defined by the species in `frames`.

Constructs an autodE `NCIComplex` from the species in `frames`
and finds the lowest energy conformation of species. Returns
a new ExtXYZ frame containing the geometry of this conformation.

The complexity of this conformer search can be modified by
changing `KineticaASE.ade.Config` before running any calculations,
using the options shown in the autodE documentation
(https://duartegroup.github.io/autodE/examples/nci.html).
"""
function autode_NCI_conformer_search(frames::Vector{Dict{String, Any}})
    @debug "Searching for conformers of reactive complex."
    mols = [frame_to_autode(frame) for frame in frames]
    sys = ade.NCIComplex(mols...)
    sys.find_lowest_energy_conformer(lmethod=ade.methods.XTB())
    n_confs_found = pylen(sys.conformers)
    @debug "$(n_confs_found) conformers found."

    sysframe = autode_to_frame(sys)
    sysframe["info"]["energy"] = pyconvert(Float64, sys.energy.to("ev").real)
    rm(joinpath(pwd(), "conformers"), recursive=true)
    return sysframe
end