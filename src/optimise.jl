"""
    geomopt!(frame::Dict{String, Any}, asecalc::Py[, fmax=0.01, maxiters=nothing])

Runs an ASE-driven geometry optimisation of the species in `frame`.

`asecalc` should be a valid ASE calculator object. This
optimisation always uses ASE's `BFGSLineSearch` optimiser,
but the maximum force `fmax` and maximum number of iterations
`maxiters` that this optimiser converges with can be
controlled with this method's keyword arguments.

Directly modifies the atomic positions and energy of the
passed in `frame`. Energies returned are in eV.
"""
function geomopt!(frame::Dict{String, Any}, asecalc::Py; fmax=0.01, maxiters=nothing)
    @debug "Starting geometry optimisation."
    atoms = frame_to_atoms(frame)
    atoms.calc = asecalc
    opt = aseopt.QuasiNewton(atoms)
    opt.run(fmax=fmax, steps=maxiters)

    frame["arrays"]["pos"] = pyconvert(Matrix, atoms.get_positions().T)
    frame["info"]["energy_ASE"] = atoms.get_potential_energy()
    frame["info"]["inertias"] = pyconvert(Vector{Float64}, atoms.get_moments_of_inertia())
    return
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
    c1 = Py(frame1["arrays"]["pos"]).to_numpy()
    c2 = Py(frame2["arrays"]["pos"]).to_numpy()
    frame1["arrays"]["pos"] = pyconvert(Matrix, rmsd.kabsch_fit(c1, c2))
    return
end