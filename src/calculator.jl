"""
"""
mutable struct ASENEBCalculator{kmType, uType, tType} <: Kinetica.AbstractKineticCalculator
    atoms_calc::Py
    neb_k
    ftol
    climb::Bool
    climb_ftol
    maxiters::Int
    interpolation::String
    n_images::Int
    remove_unconverged::Bool
    k_max::kmType
    t_unit::String
    t_mult::tType
    cached_rids::Vector{Bool}
    ts_cache::Dict{Symbol, Any}
    sd::SpeciesData
    rd::RxData
end

"""
    ASENEBCalculator(atoms_calc[, neb_k, ftol, climb, climb_ftol, maxiters, interpolation, n_images, remove_unconverged, k_max, t_unit])

Outer constructor method for ASE-driven NEB-based kinetic calculator.
"""
function ASENEBCalculator(atoms_calc::Py; neb_k=0.1, ftol=0.01, climb::Bool=true, climb_ftol=0.1,
                          maxiters=200, interpolation::String="idpp", n_images=11, remove_unconverged::Bool=true,
                          k_max::Union{Nothing, uType}=nothing, t_unit::String="s"
                         ) where {uType <: AbstractFloat}

    throw(error("atoms_calc needs to be refactored to a builder function before this can be run."))

    # Ensure autodE can get to XTB
    @assert pyconvert(Bool, ade.methods.XTB().is_available)
    
    t_mult = tconvert(t_unit, "s")
    cached_rids = Vector{Bool}()
    ts_cache = Dict{Symbol, Any}(
        :xyz => Vector{Dict{String, Any}}(),
        :reacsys_energies => Vector{Float64}(),
        :prodsys_energies => Vector{Float64}(),
        :vib_energies => Vector{Vector{Float64}}(),
        :symmetry => Vector{Int}(),
        :geometry => Vector{Int}(),
        :spin => Vector{Float64}()
    )
    sd_blank, rd_blank = init_network()
    return ASENEBCalculator(atoms_calc, neb_k, ftol, climb, climb_ftol, maxiters, 
                            interpolation, n_images, remove_unconverged, k_max, 
                            t_unit, t_mult, cached_rids, ts_cache, sd_blank, rd_blank)
end

"""
"""
function setup_network!(sd::SpeciesData{iType}, rd::RxData, calc::ASENEBCalculator) where {iType}
    # Check the cache is not out of sync.
    if count(calc.cached_rids) != length(calc.ts_cache[:atoms])
        throw(ErrorException("TS cache is out of sync! Ensure ts_cache entries match cached_rids."))
    end

    # If nothing has been cached, create an empty cache.
    if length(calc.cached_rids) == 0
        @info "No reactions cached."
        calc.cached_rids = [false for _ in 1:rd.nr]
    # If the cache is smaller than nr, extend it.
    elseif length(calc.cached_rids) < rd.nr
        @info "New reactions added to CRN, expanding reaction cache."
        calc.cached_rids = cat(calc.cached_rids, [false for _ in length(calc.cached_rids)+1:rd.nr]; dims=1)
    else
        @info "No new reactions added to CRN."
    end

    # If SpeciesData does not have caches for these properties, create them.
    if !(:vib_energies in keys(sd.cache))
        @debug "Creating empty species caches"
        sd.cache[:vib_energies] = Dict{iType, Float64}()
        sd.cache[:symmetry] = Dict{iType, Int}()
        sd.cache[:spin] = Dict{iType, Float64}()
        sd.cache[:geometry] = Dict{iType, Int}()
    end

    # Determine the species which are in reactions.
    # Generate conformers for these species, optimise and cache properties.
    active_species = unique(reduce(vcat, [rd.id_reacs; rd.id_prods]))
    @info "Searching for best conformations of all active species in CRN."
    @debug "$(active_species)/$(sd.n) species are active in the current CRN"
    for i in active_species
        if !(i in keys(sd.cache[:symmetry]))
            autode_conformer_search!(sd, i)
            geomopt!(sd.xyz[i], calc.atoms_calc)
        end
    end
    get_species_stats!(sd, refresh=true)

    # Main loop over all reactions in CRN.
    for i in 1:rd.nr
        # Skip reaction if already cached.
        if calc.cached_rids[i] continue end

        # Copy TS info if reverse reaction has been cached.
        # Reactants/products don't need to be changed because they
        # are covered by the reverse reaction.
        reverse_mapped_rxn = join(reverse(split(rd.mapped_rxns[i], ">>")), ">>")
        if reverse_mapped_rxn in rd.mapped_rxns
            reverse_idx = findfirst(==(reverse_mapped_rxn), rd.mapped_rxns)
            if calc.cached_rids[reverse_idx]
                push!(calc.ts_cache[:xyz], calc.ts_cache[:xyz][reverse_idx])
                push!(calc.ts_cache[:vib_energies], calc.ts_cache[:vib_energies][reverse_idx])
                push!(calc.ts_cache[:reacsys_energies], calc.ts_cache[:reacsys_energies][reverse_idx])
                push!(calc.ts_cache[:prodsys_energies], calc.ts_cache[:prodsys_energies][reverse_idx])
                push!(calc.ts_cache[:symmetry], calc.ts_cache[:symmetry][reverse_idx])
                push!(calc.ts_cache[:geometry], calc.ts_cache[:geometry][reverse_idx])
                push!(calc.ts_cache[:spin], calc.ts_cache[:spin][reverse_idx]) # does this actually hold?
                calc.cached_rids[i] = true
                continue
            end
        end

        # Create optimised non-covalent interacting reaction complexes for
        # reactants and products if required.
        if sum(rd.stoic_reacs[i]) > 1
            reac_frames = []
            for (j, rid) in enumerate(rd.id_reacs[i])
                for _ in 1:rd.stoic_reacs[i][j]
                    push!(reac_frames, sd.xyz[rid])
                end
            end
            reacsys = autode_NCI_conformer_search(reac_frames)
            geomopt!(reacsys, calc.atoms_calc)
        else
            reacsys = sd.xyz[rd.id_reacs[i][1]]
        end
        if sum(rd.stoic_prods[i]) > 1
            prod_frames = []
            for (j, pid) in enumerate(rd.id_prods[i])
                for _ in 1:rd.stoic_prods[i][j]
                    push!(prod_frames, sd.xyz[pid])
                end
            end
            prodsys = autode_NCI_conformer_search(prod_frames)
            geomopt!(prodsys, calc.atoms_calc)
        else
            prodsys = sd.xyz[rd.id_prods[i][1]]
        end

        # Atom map endpoints.
        reac_map, prod_map = split(rd.mapped_rxns[i], ">>")
        atom_map_frame(reac_map, reacsys)
        atom_map_frame(prod_map, prodsys)

        # Kabsch fit product system onto reactant system.
        kabsch_fit!(prodsys, reacsys)

        # Interpolate and run NEB.
        rxn_spin = determine_rxn_spin()
        images = neb(reacsys, prodsys, calc)
        # Save to caches.
        # If unconverged and removal is requested, push blank
        # entries to caches so splice! still works at the end.
        if is_neb_converged(images, calc.ftol) || !(calc.remove_unconverged)
            ts = highest_energy_frame(images)

            push!(calc.ts_cache[:xyz], ts)
            ts_sym, ts_geom = autode_frame_symmetry(ts)
            push!(calc.ts_cache[:symmetry], ts_sym)
            push!(calc.ts_cache[:geometry], ts_geom)
            rd.dH[i] = (prodsys["info"]["energy_ASE"] - reacsys["info"]["energy_ASE"]) * Constants.eV_to_kcal_per_mol
            push!(calc.ts_cache[:reacsys_energies], reacsys["info"]["energy_ASE"])
            push!(calc.ts_cache[:prodsys_energies], prodsys["info"]["energy_ASE"])
            push!(calc.ts_cache[:spin], rxn_spin)

            # Run individual vibrational analyses on reactants, products
            # and TS.
            for rid in rd.id_reacs[i]
                calc_species_vibrations!(sd, rid, calc.atoms_calc)
            end
            for pid in rd.id_prods[i]
                calc_species_vibrations!(sd, pid, calc.atoms_calc)
            end
            calc_ts_vibrations!(calc.ts_cache, i, calc.atoms_calc)
        else
            push!(calc.ts_cache[:xyz], Dict{String, Any}())
            push!(calc.ts_cache[:symmetry], -1)
            push!(calc.ts_cache[:geometry], -1)
            push!(calc.ts_cache[:reacsys_energies], 0.0)
            push!(calc.ts_cache[:prodsys_energies], 0.0)
            push!(calc.ts_cache[:spin], 0.0)
        end

        calc.cached_rids[i] = true
    end

    # Find all reactions that need to be removed from negative
    # TS symmetry entries. 
    if calc.remove_unconverged
        rem_idxs = findall(==(-1), calc.ts_cache[:symmetry])
        # Remove unconverged reactions from CRN.
        splice!(rd, calc, rem_idxs)
    end

    # Share final CRN with calculator. 
    calc.sd = sd
    calc.rd = rd

    return
end

"""
"""
function Base.splice!(calc::ASENEBCalculator, rids::Vector{Int})
    for key in keys(calc.cache)
        splice!(calc.cache[key], rids)
    end
    splice!(calc.cached_rids, rids)
end


"""
"""
function get_entropy(sd::SpeciesData, sid, T, P)
    return get_entropy(
        sd.cache[:weights][sid],
        sd.xyz[sid]["info"]["inertias"],
        sd.cache[:geometry][sid],
        sd.cache[:symmetry][sid],
        sd.cache[:spin][sid],
        sd.cache[:vib_energies][sid],
        T, P
    )
end

function get_entropy(ts_cache::Dict{Symbol, Any}, rid, mass, spin, T, P)
    return get_entropy(
        mass,
        ts_cache[:xyz][rid]["info"]["inertias"],
        ts_cache[:geometry][rid],
        ts_cache[:symmetry][rid],
        spin,
        ts_cache[:vib_energies][rid], 
        T, P
    )
end

function get_entropy(mass, inertias, geometry, symmetry, spin, vib_energies, T, P)
    S = 0.0
    
    # Translational entropy
    mass_kg = mass * ASEConstants.amu
    S_t = (2.0 * pi * mass_kg * ASEConstants.k * T / (ASEConstants.hplanck^2))^1.5
    S_t *= ASEConstants.kB * T / ASEConstants.ref_P
    S_t = ASEConstants.kB * (log10(S_t) + 2.5)
    S += S_t

    # Rotational entropy
    if geometry == 1
        inertias_conv = inertias * ASEConstants.amu / (10.0^10)^2
        inertia = maximum(inertias_conv)
        S_r = 8.0 * pi^2 * inertia * ASEConstants.k * T / symmetry / ASEConstants.hplanck^2
        S_r = ASEConstants.kB * (log10(S_r) + 1.0)
    elseif geometry == 2
        inertias_conv = inertias * ASEConstants.amu / (10.0^10)^2
        S_r = sqrt(pi * prod(inertias_conv)) / symmetry
        S_r *= (8.0 * pi^2 * ASEConstants.k * T / ASEConstants.hplanck^2)^1.5
        S_r = ASEConstants.kB * (log10(S_r) + 1.5)
    else
        S_r = 0.0
    end
    S += S_r

    # Electronic entropy
    S += ASEConstants.kB * log10(2*spin + 1)

    # Vibrational entropy
    kT = ASEConstants.kB * T
    S_v = 0.0
    for e in vib_energies
        x = e/kT
        S_v += x / (exp(x) - 1.0) - log10(1.0 - exp(-x))
    end
    S += S_v * ASEConstants.kB

    # Pressure correction to translational entropy
    S_p = -ASEConstants.kB * log10(P / ASEConstants.ref_P)
    S += S_p

    return S
end

"""
"""
function get_enthalpy(sd::SpeciesData, sid, T)
    return get_enthalpy(sd.xyz[sid]["info"]["energy_ASE"], sd.cache[:vib_energies][sid], sd.cache[:geometry][sid], T)
end

function get_enthalpy(ts_cache::Dict{Symbol, Any}, rid, T)
    return get_enthalpy(ts_cache[:xyz][rid]["info"]["energy_ASE"], ts_cache[:vib_energies][rid], ts_cache[:geometry][rid], T)
end

function get_enthalpy(energy, vib_energies, geometry, T)
    H = 0.0
    H += energy

    # Add ZPE correction.
    for e in vib_energies
        H += 0.5 * e
    end

    # Translational heat capacity
    H += 1.5 * ASEConstants.kB * T

    # Rotational heat capacity
    if geometry == 1
        H += ASEConstants.kB * T
    elseif geometry == 2
        H += ASEConstants.kB * 1.5 * T
    end

    # Vibrational heat capacity
    kT = ASEConstants.kB * T
    for e in vib_energies
        H += e / (exp(e/kT) - 1)
    end

    H += ASEConstants.kB * T
    return H
end

# Dispatched with k_max awareness.
"""
    calculator(; T, P)
"""
function (calc::ASENEBCalculator{uType, uType, tType})(; T::Number, P::Number) where {uType, tType}
    throw(ErrorException(""))
end

# Dispatched without k_max awareness.
function (calc::ASENEBCalculator{Nothing, uType, tType})(; T::Number, P::Number) where {uType, tType}
    dS = zeros(calc.rd.nr)
    dH = zeros(calc.rd.nr)
    for rid in 1:calc.rd.nr
        S_reacs = 0.0
        H_reacs = 0.0
        mass_ts = 0.0
        for (i, sid) in enumerate(calc.rd.id_reacs[rid])
            spec_stoic = calc.rd.stoic_reacs[rid][i]
            mass_ts += spec_stoic * calc.sd.cache[:weights][sid]
            S_reacs += spec_stoic * get_entropy(calc.sd, sid, T, P)
            H_reacs += spec_stoic * get_enthalpy(calc.sd, sid, T)
        end


    end

    k = Constants.k_b*T/Constants.h .* exp.(dS/Constants.R) .* exp.(-dH/(Constants.R*T))
    return k
end

function has_conditions(::ASENEBCalculator, symbols::Vector{Symbol})
    return all([sym in [:T, :P] for sym in symbols])
end

allows_continuous(::ASENEBCalculator) = false