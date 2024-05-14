"""
"""
function get_rxn_mult(reacsys, prodsys)
    n_reacs = reacsys["info"]["n_species"]
    reac_mult = reacsys["info"]["mult"]
    n_prods = prodsys["info"]["n_species"]
    prod_mult = prodsys["info"]["mult"]

    # Given an addition (2->1) or a dissociation (1->2), pick the lower mult.
    if n_reacs > n_prods
        return prod_mult
    elseif n_reacs < n_prods
        return reac_mult
    # Given a substitution (2->2) or a rearrangement (1->1), 
    # the mult should always be balanced.
    else
        if reac_mult != prod_mult
            throw(ErrorException("Reaction has equal n_reacs and n_prods but mults do not match."))
        else
            return reac_mult
        end
    end
end


"""
    neb(reacsys, prodsys, calc::ASENEBCalculator[, calcdir="./", kwargs...])

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
function neb(reacsys, prodsys, calc::ASENEBCalculator; calcdir="./", kwargs...)
    @info "Running $(calc.climb ? "CI-" : "")NEB calculation"
    images = [
        [frame_to_atoms(reacsys, reacsys["info"]["formal_charges"]) for _ in 1:calc.n_images]; 
        [frame_to_atoms(prodsys, prodsys["info"]["formal_charges"])]
    ]

    rmult = get_rxn_mult(reacsys, prodsys)
    if calc.parallel
        for image in images
            image.calc = calc.calc_builder(calcdir, rmult, reacsys["info"]["chg"], kwargs...)
        end
    else
        shared_calc = calc.calc_builder(calcdir, rmult, reacsys["info"]["chg"], kwargs...)
        for image in images
            image.calc = shared_calc
        end
    end
    images = pylist(images)
    neb = aseneb.NEB(images, k=calc.neb_k, parallel=calc.parallel, allow_shared_calculator=(!calc.parallel)) 

    @debug "Interpolating reaction path with method: $(calc.interpolation)"
    if calc.interpolation in ["linear", "idpp"]
        neb.interpolate(method=calc.interpolation)
    else
        throw(ErrorException("Unknown interpolation method. must be one of [\"linear\", \"idpp\"]"))
    end
    aseio.write(joinpath(calcdir, "interp.traj"), images)

    # opt = aseopt.FIRE(neb)
    opt = aseneb.NEBOptimizer(neb, verbose=1)
    if calc.climb
        @debug "Running NEB to tolerance of $(calc.climb_ftol) before enabling CI"
        conv = opt.run(fmax=calc.climb_ftol, steps=calc.maxiters)
        conv = pyconvert(Bool, pybuiltins.bool(conv))
        if conv
            @debug "Running CI-NEB to tolerance of $(calc.ftol)"
            neb.climb = true
            conv = opt.run(fmax=calc.ftol, steps=calc.maxiters)
            conv = pyconvert(Bool, pybuiltins.bool(conv))
        end
    else
        @debug "Running NEB to tolerance of $(calc.ftol)"
        conv = opt.run(fmax=calc.ftol, steps=calc.maxiters)
        conv = pyconvert(Bool, pybuiltins.bool(conv))
    end
    aseio.write(joinpath(calcdir, "neb_final.traj"), images)

    final_fmax = pyconvert(Float64, opt.get_residual())
    if conv
        @info "NEB converged (fmax = $(final_fmax))"
    else
        @info "NEB not converged (fmax = $(final_fmax))"
    end

    return images, conv
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
    ts["info"]["formal_charges"] = pyconvert(Vector{Float64}, images[ts_idx-1].get_initial_charges())
    return ts
end