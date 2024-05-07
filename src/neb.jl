function determine_rxn_spin()
    throw(error("Not yet implemented!"))
end


"""
    neb(reacsys, prodsys, calc::ASENEBCalculator)

Interpolates and runs (CI)NEB for the reaction defined by endpoints `reacsys` and `prodsys`.

Given ExtXYZ frames for reactants and products `reacsys`
and `prodsys` respectively, performs an interpolation
to generate an initial reaction path and then optimises
it with ASE's NEB implementation.

Interpolation scheme is selected by `calc.interpolation`, 
where `"linear"` and `"idpp"` are the currently available
options.

NEB optimisation is performed by the FIRE optimiser. If
`calc.climb=false`, runs a regular NEB calculation until
maximum force is less than `calc.ftol`. If `calc.climb=true`,
runs a regular NEB calculation until maximum force is less
than `calc.climb_ftol`, then enables CINEB and runs until
maximum force is less than `calc.ftol`. Each optimisation
only runs until `calc.maxiters` iterations have elapsed.

Returns the optimised NEB path as a Python list of ASE
Atoms objects.
"""
function neb(reacsys, prodsys, calc::ASENEBCalculator)
    @debug "Running $(calc.climb ? "CI-" : "")NEB calculation"
    images = [[frame_to_atoms(reacsys) for _ in 1:calc.n_images]; [frame_to_atoms(prodsys)]]
    for image in images
        image.calc = calc.atoms_calc
    end
    images = pylist(images)

    @debug "Interpolating reaction path with method: $(calc.interpolation)"
    if calc.interpolation in ["linear", "idpp"]
        neb = aseneb.NEB(images, k=calc.neb_k, allow_shared_calculator=true) 
        neb.interpolate(method=calc.interpolation)
    else
        throw(ErrorException("Unknown interpolation method. must be one of [\"linear\", \"idpp\"]"))
    end

    opt = aseopt.FIRE(neb)
    if calc.climb
        @debug "Running NEB to tolerance of $(calc.climb_ftol) before enabling CI"
        opt.run(fmax=calc.climb_ftol, steps=calc.maxiters)
        @debug "Running CI-NEB to tolerance of $(calc.ftol)"
        neb.climb = true
        opt.run(fmax=calc.ftol, steps=calc.maxiters)
    else
        @debug "Running NEB to tolerance of $(calc.ftol)"
        opt.run(fmax=calc.ftol, steps=calc.maxiters)
    end

    return images
end


"""
    is_neb_converged(images::Py, ftol)

Checks whether a NEB calculation converged to the requested force tolerance `ftol`.

Mostly exists to output a debug message.
"""
function is_neb_converged(images::Py, ftol)
    final_fmax = pyconvert(Float64, aseneb.NEBTools(images).get_fmax())
    if final_fmax < ftol
        @debug "NEB converged (fmax = $(final_fmax))"
        return true
    else
        @debug "NEB not converged (fmax = $(final_fmax))"
        return false
    end
end

"""
    highest_energy_frame(images::Py)

Finds the highest energy NEB image in `images`, returns as a frame.

`images` should be a Python list of ASE Atoms objects
representing a NEB path (see the `neb(reacsys, prodsys, atoms_calc)`
method). Locates the highest energy image, converts it
to an ExtXYZ frame and returns this frame.
"""
function highest_energy_frame(images::Py)
    energies = [pyconvert(Float64, image.get_potential_energy()) for image in images]
    ts_idx = argmax(energies)
    @debug "TS found at image $(ts_idx)/$(pylen(images))"
    inertias = pyconvert(Vector{Float64}, images[ts_idx-1].get_moments_of_inertia())
    ts = atoms_to_frame(images[ts_idx-1], energies[ts_idx], inertias)
    return ts
end