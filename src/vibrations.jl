"""
    calc_species_vibrations(sd::SpeciesData, sid, atoms_calc::Py[, refresh=false])

Calculates vibrational energies of a species.

Given a species with an optimised geometry in `sd.xyz[sid]`,
uses ASE to calculate vibrational energies with a finite
difference approximation of its Hessian, using the ASE
calculator `atoms_calc`.

Writes the vibrational energies in eV as an array into 
`sd.cache[:vib_energies]`. Defaults to not overwriting 
already calculated vibrational energies for a species,
instead returnong immediately. If recalculations are
desired, use `refresh=true`. 
"""
function calc_species_vibrations!(sd::SpeciesData, sid, atoms_calc::Py; refresh=false)
    if sid in keys(sd.cache[:vib_energies]) && !(refresh)
        @debug "Species $sid has vibrations cached, skipping."
        return
    end

    if sd.cache[:geometry][sid] == 0
        @debug "Species $sid is monoatomic, skipping vibrational analysis."
        sd.cache[:vib_energies][sid] = []
    end

    atoms = frame_to_atoms(sd.xyz[sid])
    atoms.calc = atoms_calc

    vibdir = joinpath(pwd(), "vib")
    if !(isdir(vibdir)) mkdir(vibdir) end
    vib = asevib.Vibrations(atoms)
    vib.run()
    vib_energies = vib.get_energies()
    if sd.cache[:geometry][sid] == 1
        vib_energies = vib_energies[pyslice(-(3*sd.xyz[sid]["N_atoms"] - 6), pylen(vib_energies))]
    elseif sd.cache[:geometry][sid] = 2
        vib_energies = vib_energies[pyslice(-(3*sd.xyz[sid]["N_atoms"] - 5), pylen(vib_energies))]
    else
        throw(ErrorException("Unknown geometry for species $sid: $(sd.cache[:geometry][sid])"))
    end
    jlve = pyconvert(Vector{ComplexF64}, vib_energies)
    if any([z.im > 1e-8 for z in jlve])
        throw(ErrorException("Imaginary frequency detected in geometry of species $sid."))
    else
        jlve = [z.re for z in jlve]
    end
    sd.cache[:vib_energies][sid] = jlve
    rm(vibdir, recursive=true)
    return
end

"""
    calc_ts_vibrations(ts_cache::Dict{Symbol, Any}, rid, atoms_calc::Py)

Calculates vibrational energies of a transition state.

Given a TS with an optimised geometry in `ts_cache[:xyz][rid]`,
uses ASE to calculate vibrational energies with a finite 
difference apporoximation of its Hessian, using the ASE
calculator `atoms_calc`.

Writes the vibrational energies in eV as an array into 
`ts_cache[:vib_energies]`.
"""
function calc_ts_vibrations!(ts_cache::Dict{Symbol, Any}, rid, atoms_calc::Py)
    atoms = frame_to_atoms(ts_cache[:xyz][rid])
    atoms.calc = atoms_calc

    vibdir = joinpath(pwd(), "vib")
    if !(isdir(vibdir)) mkdir(vibdir) end
    vib = asevib.Vibrations(atoms)
    vib.run()
    vib_energies = vib.get_energies()
    if ts_cache[:geometry][rid] == 1
        vib_energies = vib_energies[pyslice(-(3*ts_cache[:xyz][rid]["N_atoms"] - 5), pylen(vib_energies))]
    elseif ts_cache[:geometry][rid] == 2
        vib_energies = vib_energies[pyslice(-(3*ts_cache[:xyz][rid]["N_atoms"] - 6), pylen(vib_energies))]
    else
        throw(ErrorException("Unknown geometry for TS $rid: $(ts_cache[:geometry][rid])"))
    end
    jlve = pyconvert(Vector{ComplexF64}, vib_energies)
    if any([z.im > 1e-8 for z in jlve])
        throw(ErrorException("Imaginary frequency detected in goemetry of TS $rid."))
    else
        jlve = [z.re for z in jlve]
    end
    push!(ts_cache[:vib_energies], jlve)
    rm(vibdir, recursive=true)
    return
end

