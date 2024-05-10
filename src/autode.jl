
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
            pyconvert(Float64, ademol.atoms[i].coord.x),
            pyconvert(Float64, ademol.atoms[i].coord.y),
            pyconvert(Float64, ademol.atoms[i].coord.z)
        ] for i in 0:na-1
    ])
    
    info = isnothing(info_dict) ? Dict{String, Any}() : info_dict
    arrays = Dict{String, Any}("species" => symbols, "pos" => coords)
    frame = Dict{String, Any}("N_atoms" => na, "info" => info, "arrays" => arrays)
    return frame
end

"""
    frame_to_autode(frame::Dict{String, Any}[, mult::Int=1, chg::Int=0])

Converts an ExtXYZ frame to an autodE `Molecule`.

Writes a temporary xyz file from `frame` and reads it back in
with autodE. There doesn't seem to be a way around writing to
disk here, as autodE only detects an xyz input when reading a
file ending in '.xyz', which can't be replicated with an in-
memory IO buffer.

Optionally allows specification of spin multiplicity and charge
through the `mult` and `chg` keyword arguments respectively.
"""
function frame_to_autode(frame::Dict{String, Any}; mult::Int=1, chg::Int=0)
    f = joinpath(tempdir(), "frame.xyz")
    write_frame(f, frame)
    mol = ade.Molecule(f, charge=chg, mult=mult)
    rm(f)
    return mol
end

"""
    autode_conformer_search!(sd::SpeciesData, i)

Performs a conformer search for the species at `sd.xyz[i]`, updating its geometry.

Constructs an initial guess of species geometry from its SMILES,
then runs an autodE conformer search with xTB as the energetic
driver. Finds the lowest energy conformer and writes it back
to `sd.xyz[i]`.

Requires spin multiplicity and charge for the given species to be
cached in `sd.cache[:mult][i]` and `sd.cache[:charge][i]` 
respectively, which can be achieved by calling `get_mult!` and 
`get_charge!`. Failure to do so will result in an error.

Also populates `sd.cache` with autodE-derived values for symmetry
number and geometry for later use in TST calculations.
"""
function autode_conformer_search!(sd::SpeciesData, i)
    if !(i in keys(sd.cache[:mult])) || !(i in keys(sd.cache[:charge]))
        throw(KeyError("Missing multiplicity and/or charge in cache for SID $i."))
    end

    mol = ade.Molecule(smiles=sd.toStr[i], mult=sd.cache[:mult][i], charge=sd.cache[:charge][i])
    if sd.xyz[i]["N_atoms"] > 2
        @debug "Searching for conformers of species $i: $(sd.toStr[i]) (mult = $(sd.cache[:mult][i]), charge = $(sd.cache[:charge][i]))"
        mol.find_lowest_energy_conformer()
        n_confs_found = pylen(mol.conformers)
        @debug "$(n_confs_found) conformers found."
    else
        @debug "Species $i ($(sd.toStr[i])) too small for conformer search."
        mol.optimise(method=ade.methods.XTB())
    end

    @debug "Writing lowest energy conformer to SpeciesData."
    sd.xyz[i] = autode_to_frame(mol; info_dict=sd.xyz[i]["info"])
    sd.xyz[i]["info"]["energy"] = pyconvert(Float64, mol.energy.to("ev").real)
    
    sd.cache[:symmetry][i] = pyconvert(Int, mol.symmetry_number)
    if pyconvert(Int, mol.n_atoms) == 1
        sd.cache[:geometry][i] = 0
    else
        if pyconvert(Bool, mol.is_linear())
            sd.cache[:geometry][i] = 1
        else
            sd.cache[:geometry][i] = 2
        end
    end

    # rm(joinpath(pwd(), "conformers"), recursive=true)
end

"""
    autode_NCI_conformer_search(sd::SpeciesData, sids[, name])
    
Performs a search for the lowest energy non-covalent interacting reaction complex defined by the species in `sd` at species IDs `sids`.

Constructs an autodE `NCIComplex` from the provided species
and finds the lowest energy conformation of species. Returns
a new ExtXYZ frame containing the geometry of this conformation.

This frame is additionally tagged with information about the
complex, such as its spin multiplicity, charge and xTB energy.
This information can be found in the resulting `frame`'s "info"
dictionary.

This method should only be performed following calls to 
`autode_conformer_search!`, as it both generates correct species
geometries and also populates `sd.cache` with neccessary
information about the spins and charges of species.

The complexity of this conformer search can be modified by
changing `KineticaASE.ade.Config` before running any calculations,
using the options shown in the autodE documentation
(https://duartegroup.github.io/autodE/examples/nci.html).
"""
function autode_NCI_conformer_search(sd::SpeciesData, sids; name="complex")
    @debug "Searching for conformers of reactive complex."
    mols = [frame_to_autode(sd.xyz[i], mult=sd.cache[:mult][i], chg=sd.cache[:charge][i]) for i in sids]
    sys = ade.NCIComplex(mols..., name=name, do_init_translation=true)

    # smisys = ade.NCIComplex([ade.Molecule(smiles=sd.toStr[i], mult=sd.cache[:mult][i], chg=sd.cache[:charge][i]) for i in sids]...)
    # ade.mol_graphs.make_graph(smisys)

    sys._generate_conformers()
    sys.conformers.optimise(method=ade.methods.XTB())
    sys.conformers.remove_no_energy()
    sys.conformers.prune_diff_graph(graph=sys.graph)
    if pylen(sys.conformers) < 1
        throw(ErrorException("All generated conformers break molecular graph."))
    end
    sys.conformers.prune()
    sys._set_lowest_energy_conformer()
    # sys.find_lowest_energy_conformer(lmethod=ade.methods.XTB())
    n_confs_found = pylen(sys.conformers)
    @debug "$(n_confs_found) conformers found."

    sysframe = autode_to_frame(sys; info_dict=Dict{String, Any}(
        "mult" => pyconvert(Int, sys.mult),
        "chg" => pyconvert(Int, sys.charge),
        "energy" => pyconvert(Float64, sys.energy.to("ev").real),
        "n_species" => length(sids)
    ))
    for f in glob("./*_conf*")
        rm(f)
    end
    return sysframe
end


"""
"""
function autode_frame_symmetry(frame::Dict{String, Any})
    mol = frame_to_autode(frame)
    sym = pyconvert(Int, mol.symmetry_number)
    if pyconvert(Bool, mol.is_linear())
        geom = 1
    else
        geom = 2
    end
    return sym, geom
end