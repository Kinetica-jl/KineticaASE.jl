"""
"""
function get_mult(sd, i)
    smi = sd.toStr[i]
    mol = rdChem.MolFromSmiles(smi)
    n_radical_electrons = 0
    for atom in mol.GetAtoms()
        n_radical_electrons += pyconvert(Int, atom.GetNumRadicalElectrons())
    end
    mult = n_radical_electrons + 1
    return mult
end

"""
"""
get_mult!(sd, i) = sd.cache[:mult][i] = get_mult(sd, i)
    

"""
"""
get_charge(sd, i) = pyconvert(Int, pybel.readstring("can", sd.toStr[i]).charge)
get_charge!(sd, i) = sd.cache[:charge][i] = get_charge(sd, i)


"""
"""
function get_formal_charges(amsmi::String)
    mol = pybel.readstring("smi", amsmi)
    formal_charges = [pyconvert(Int, atom.formalcharge) for atom in mol.atoms]
    return formal_charges
end
function get_formal_charges(sd::SpeciesData, i)
    amsmi = atom_map_smiles(sd.xyz[i], sd.toStr[i])
    return get_formal_charges(amsmi)
end
get_formal_charges!(sd::SpeciesData, i) = sd.cache[:formal_charges][i] = get_formal_charges(sd, i)

"""
    geomopt!(frame::Dict{String, Any}, calc_builder[, mult::Int=1, chg::Int=0, formal_charges=nothing, fmax=0.01, maxiters=nothing, kwargs...])
    geomopt!(sd::SpeciesData, i, calc_builder[, fmax=0.01, maxiters=nothing, kwargs...])

Runs an ASE-driven geometry optimisation of the species in `frame`.

Can be run directly from a `frame`, or a `frame` can be extracted
from `sd.xyz[i]` in theh case of the second method. With this method,
formal charges, total charge and spin multiplicity are assumed
to have been calculated and cached in `sd.cache`. When running
directly from a `frame`, this information must be passed manually.

`calc_builder` should be a struct with a functor that returns
a correctly constructed ASE calculator for the system at hand.
While ASE can handle many system-specific calculator details
from an `Atoms` object, quantities such as spin multiplicity
and sometimes charge must be input separately. For this reason,
the `calc_builder` functor must take `mult::Int` and `charge::Int`
as its first two arguments. Any other arguments can be passed
via this method's `kwargs`.

Some ASE calculators handle charged species at the `Atoms` level.
These require an array of formal charges on each atom to be
given during `Atoms` construction. If `formal_charges` is
provided, this will occur. If not, all formal charges will be 
assumed to be zero. 

IMPORTANTLY, any `formal_charges` array MUST match the atom
ordering of the provided `frame`. This can by calling 
`get_formal_charges` on an atom-mapped SMILES from `atom_map_smiles`.

This optimisation always uses ASE's `BFGSLineSearch` optimiser,
but the maximum force `fmax` and maximum number of iterations
`maxiters` that this optimiser converges with can be
controlled with this method's keyword arguments.

Directly modifies the atomic positions and energy of the
passed in `frame`. Energies returned are in eV. Returns a
boolean for whether the optimisation was a conv.
"""
function geomopt!(sd::SpeciesData, i, calc_builder; calcdir::String="./", fmax=0.01, maxiters=nothing, kwargs...)
    frame = sd.xyz[i]
    conv = geomopt!(frame, calc_builder; calcdir=calcdir, mult=sd.cache[:mult][i], chg=sd.cache[:charge][i],
                       formal_charges=sd.cache[:formal_charges][i], fmax=fmax, maxiters=maxiters,
                       kwargs...)
    sd.xyz[i] = frame
    return conv
end

function geomopt!(frame::Dict{String, Any}, calc_builder; 
                  calcdir::String="./", mult::Int=1, chg::Int=0, 
                  formal_charges=nothing, fmax=0.01, maxiters=nothing,
                  kwargs...)
    @debug "Starting geometry optimisation."
    atoms = frame_to_atoms(frame, formal_charges)
    atoms.calc = calc_builder(calcdir, mult, chg, kwargs...)
    init_energy = pyconvert(Float64, atoms.get_potential_energy())
    init_inertias = pyconvert(Vector{Float64}, atoms.get_moments_of_inertia())

    opt = aseopt.QuasiNewton(atoms); conv = false
    try
        conv = opt.run(fmax=fmax, steps=maxiters)
        conv = pyconvert(Bool, pybuiltins.bool(conv))
    catch err
        conv = false
    end

    if conv
        @debug "Geometry optimisation complete."
        frame["arrays"]["pos"] = pyconvert(Matrix, atoms.get_positions().T)
        frame["info"]["energy_ASE"] = pyconvert(Float64, atoms.get_potential_energy())
        frame["info"]["inertias"] = pyconvert(Vector{Float64}, atoms.get_moments_of_inertia())
    else
        @debug "Geometry optimisation failed."
        frame["info"]["energy_ASE"] = init_energy
        frame["info"]["inertias"] = init_inertias
    end
    return conv
end


"""
"""
function safe_geomopt!(frame::Dict{String, Any}, calc_builder; 
                       calcdir::String="./", mult::Int=1, chg::Int=0, 
                       formal_charges=nothing, fmax=0.01, maxiters=nothing,
                       kwargs...)
    @debug "Starting geometry optimisation (safe)."
    ademol_orig = frame_to_autode(frame; mult=mult, chg=chg)

    atoms = frame_to_atoms(frame, formal_charges)
    atoms.calc = calc_builder(calcdir, mult, chg, kwargs...)
    init_energy = pyconvert(Float64, atoms.get_potential_energy())
    init_inertias = pyconvert(Vector{Float64}, atoms.get_moments_of_inertia())

    opt = aseopt.QuasiNewton(atoms); conv = false
    try
        conv = opt.run(fmax=fmax, steps=maxiters)
        conv = pyconvert(Bool, pybuiltins.bool(conv))
    catch err
        conv = false
    end

    if conv
        @debug "Geometry optimisation complete."
        optframe = atoms_to_frame(atoms, pyconvert(Float64, atoms.get_potential_energy()), 
                                pyconvert(Vector{Float64}, atoms.get_moments_of_inertia()))
        ademol_opt = frame_to_autode(optframe; mult=mult, chg=chg)
        if !pyconvert(Bool, ade.mol_graphs.is_isomorphic(ademol_orig.graph, ademol_opt.graph))
            @warn "Optimised geometry invalidates molecular graph, reverting to unoptimised geometry."
            frame["info"]["energy_ASE"] = init_energy
            frame["info"]["inertias"] = init_inertias
            conv = false
        else
            frame["arrays"]["pos"] = pyconvert(Matrix, atoms.get_positions().T)
            frame["info"]["energy_ASE"] = pyconvert(Float64, atoms.get_potential_energy())
            frame["info"]["inertias"] = pyconvert(Vector{Float64}, atoms.get_moments_of_inertia())
        end
    else
        @warn "Geometry optimisation failed."
        frame["info"]["energy_ASE"] = init_energy
        frame["info"]["inertias"] = init_inertias
    end
    return conv
end

"""
    kabsch_fit!(frame1::Dict{String, Any}, frame2::Dict{String, Any})

Computes the maximum overlap of atoms within `frame1` to `frame2`.

Modifies the atomic positions in `frame1` to be as close as
possible to those in `frame2` through a combination of
translation and rotation. Uses the Kabsch algorithm, as
implemented in the Python package 'rmsd'.
"""
function kabsch_fit!(frame1::Dict{String, Any}, frame2::Dict{String, Any})
    c1 = Py(frame1["arrays"]["pos"]).to_numpy().T
    c2 = Py(frame2["arrays"]["pos"]).to_numpy().T
    frame1["arrays"]["pos"] = pyconvert(Matrix, rmsd.kabsch_fit(c1, c2).T)
    return
end


"""
"""
function get_hydrogen_idxs(amsmi::String)
    at_end = false
    i = 1
    nc = length(amsmi)
    hidxs = [Int[]]
    while !at_end
        if amsmi[i] == '['
            sym = amsmi[i+1]
            idx = parse(Int, string(amsmi[i+3]))
            if sym == 'H'
                push!(hidxs[end], idx)
            end
            i += 5
        elseif amsmi[i] == '.'
            push!(hidxs, Int[])
            i += 1
        else
            i += 1
        end
        if i > nc at_end = true end
    end
    return hidxs
end


"""
"""
function permute_hydrogens!(frame1::Dict{String, Any}, hidxs::Vector{Vector{Int}}, frame2::Dict{String, Any})
    c1 = Py(frame1["arrays"]["pos"]).to_numpy().T
    c2 = Py(frame2["arrays"]["pos"]).to_numpy().T

    if length(reduce(vcat, hidxs)) > 1
        best_pos = c1.copy()
        best_rmsd = pyconvert(Float64, rmsd.kabsch_rmsd(best_pos, c2))
        swapping = true
        while swapping
            has_swapped = false
            for hidxs_mol in hidxs
                if length(hidxs_mol) < 2 continue end
                for i in 1:length(hidxs_mol)-1
                    for j in 2:length(hidxs_mol)
                        swap_pos = best_pos.copy()
                        swap_pos[hidxs_mol[i]-1] = best_pos[hidxs_mol[j]-1]
                        swap_pos[hidxs_mol[j]-1] = best_pos[hidxs_mol[i]-1]
                        swap_rmsd = pyconvert(Float64, rmsd.kabsch_rmsd(swap_pos, c2))
                        if swap_rmsd < best_rmsd
                            @debug "Swapped H$(hidxs_mol[i]) for H$(hidxs_mol[j])"
                            best_pos = swap_pos
                            best_rmsd = swap_rmsd
                            has_swapped = true
                        end
                    end
                end
            end
            if !has_swapped
                swapping = false
            end
        end
        c1 = rmsd.kabsch_fit(best_pos, c2)
    end

    frame1["arrays"]["pos"] = pyconvert(Matrix, c1.T)
    return
end