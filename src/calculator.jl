"""
"""
mutable struct ASENEBCalculator{kmType, tType} <: Kinetica.AbstractKineticCalculator
    calc_builder
    calcdir_head::String
    neb_k
    ftol
    climb::Bool
    climb_ftol
    maxiters::Int
    interpolation::String
    n_images::Int
    parallel::Bool
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
    ASENEBCalculator(calc_builder, calcdir_head[, neb_k, ftol, climb, climb_ftol, maxiters, interpolation, n_images, parallel, remove_unconverged, k_max, t_unit])

Outer constructor method for ASE-driven NEB-based kinetic calculator.
"""
function ASENEBCalculator(calc_builder, calcdir_head; neb_k=0.1, ftol=0.01, climb::Bool=true, climb_ftol=0.1,
                          maxiters=500, interpolation::String="idpp", n_images=11, parallel::Bool=false, 
                          remove_unconverged::Bool=true, k_max::Union{Nothing, uType}=nothing, 
                          t_unit::String="s") where {uType <: AbstractFloat}

    # Ensure autodE can get to XTB
    @assert pyconvert(Bool, ade.methods.XTB().is_available)

    # Ensure calc_builder doesn't error for basic construction
    test_calc = calc_builder(mktempdir(), 1, 0)
    @assert typeof(test_calc) == Py

    # Create main calculation directory
    calcdir_head = abspath(calcdir_head)
    if !isdir(calcdir_head) mkpath(calcdir_head) end
    
    t_mult = tconvert(t_unit, "s")
    cached_rids = Vector{Bool}()
    ts_cache = Dict{Symbol, Any}(
        :xyz => Vector{Dict{String, Any}}(),
        :reacsys_energies => Vector{Float64}(),
        :prodsys_energies => Vector{Float64}(),
        :vib_energies => Vector{Vector{Float64}}(),
        :symmetry => Vector{Int}(),
        :geometry => Vector{Int}(),
        :mult => Vector{Int}(),
        :charge => Vector{Int}()
    )
    sd_blank, rd_blank = init_network()
    return ASENEBCalculator(calc_builder, calcdir_head, neb_k, ftol, climb, climb_ftol, maxiters, 
                            interpolation, n_images, parallel, remove_unconverged, k_max, 
                            t_unit, t_mult, cached_rids, ts_cache, sd_blank, rd_blank)
end

"""
"""
function Kinetica.setup_network!(sd::SpeciesData{iType}, rd::RxData, calc::ASENEBCalculator) where {iType}
    # Check the cache is not out of sync.
    if count(calc.cached_rids) != length(calc.ts_cache[:xyz])
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
        sd.cache[:vib_energies] = Dict{iType, Vector{Float64}}()
        sd.cache[:symmetry] = Dict{iType, Int}()
        sd.cache[:mult] = Dict{iType, Int}()
        sd.cache[:charge] = Dict{iType, Int}()
        sd.cache[:formal_charges] = Dict{iType, Vector{Int}}()
        sd.cache[:geometry] = Dict{iType, Int}()
    end

    # Determine the species which are in reactions.
    # Generate conformers for these species, optimise and cache properties.
    currdir = pwd()
    specoptdir_head = joinpath(calc.calcdir_head, "species_opts")
    if !isdir(specoptdir_head) mkdir(specoptdir_head) end
    active_species = unique(reduce(vcat, [rd.id_reacs; rd.id_prods]))
    @info "Searching for best conformations of all active species in CRN."
    @debug "$(length(active_species))/$(sd.n) species are active in the current CRN"
    for i in active_species
        if !("energy_ASE" in keys(sd.xyz[i]["info"]))
            specoptdir = joinpath(specoptdir_head, "spec_$(lpad(i, 6, "0"))")
            if !isdir(specoptdir) mkdir(specoptdir) end; cd(specoptdir)
            get_mult!(sd, i)
            get_charge!(sd, i) 
            autode_conformer_search!(sd, i)
            get_formal_charges!(sd, i)
            conv = geomopt!(sd, i, calc.calc_builder; maxiters=calc.maxiters)
            if !conv
                @warn "Optimisation of species $i ($(sd.toStr[i])) failed to converge!"
            end
            cd(currdir)
        end
    end
    get_species_stats!(sd, refresh=true)

    # Main loop over all reactions in CRN.
    nebdir_head = joinpath(calc.calcdir_head, "nebs")
    if !isdir(nebdir_head) mkdir(nebdir_head) end
    for i in 1:rd.nr
        # Skip reaction if already cached.
        if calc.cached_rids[i] continue end
        @info "---------------------------------\nReaction $i\n---------------------------------"
        @info format_rxn(sd, rd, i)

        # Copy TS info if reverse reaction has been cached.
        # Reactants/products don't need to be changed because they
        # are covered by the reverse reaction.
        reverse_rhash = get_reverse_rhash(sd, rd, i)
        if reverse_rhash in rd.rhash
            reverse_idx = findfirst(==(reverse_rhash), rd.rhash)
            if calc.cached_rids[reverse_idx] && calc.ts_cache[:symmetry][reverse_idx] > -1 
                push!(calc.ts_cache[:xyz], calc.ts_cache[:xyz][reverse_idx])
                push!(calc.ts_cache[:vib_energies], calc.ts_cache[:vib_energies][reverse_idx])
                push!(calc.ts_cache[:reacsys_energies], calc.ts_cache[:reacsys_energies][reverse_idx])
                push!(calc.ts_cache[:prodsys_energies], calc.ts_cache[:prodsys_energies][reverse_idx])
                push!(calc.ts_cache[:symmetry], calc.ts_cache[:symmetry][reverse_idx])
                push!(calc.ts_cache[:geometry], calc.ts_cache[:geometry][reverse_idx])
                push!(calc.ts_cache[:mult], calc.ts_cache[:mult][reverse_idx]) # does this actually hold?
                push!(calc.ts_cache[:charge], calc.ts_cache[:charge][reverse_idx])
                @info "Found reverse reaction in cache, skipping."
                @info ""
                calc.cached_rids[i] = true
                continue
            end
        end

        # Create a unique directory that is independent of reaction numbering.
        nebdir = joinpath(nebdir_head, bytes2hex(rd.rhash[i]))
        if !isdir(nebdir) mkdir(nebdir) end; cd(nebdir)
        @info "Running calculations in $nebdir"

        # Create optimised non-covalent interacting reaction complexes for
        # reactants and products if required.
        if sum(rd.stoic_reacs[i]) > 1
            reac_sids = []
            for (j, sid) in enumerate(rd.id_reacs[i])
                for _ in 1:rd.stoic_reacs[i][j]
                    push!(reac_sids, sid)
                end
            end
            reacsys = autode_NCI_conformer_search(sd, reac_sids; name="reacsys")
            reacsys["info"]["n_species"] = length(reac_sids)
            reacsys_smi = join([sd.toStr[sid] for sid in reac_sids], ".")
            formal_charges = get_formal_charges(atom_map_smiles(reacsys, reacsys_smi))
            safe_geomopt!(reacsys, calc.calc_builder; calcdir=nebdir, mult=reacsys["info"]["mult"], 
                          chg=reacsys["info"]["chg"], formal_charges=formal_charges, maxiters=calc.maxiters)
        else
            sid = rd.id_reacs[i][1]
            reacsys = sd.xyz[sid]
            reacsys["info"]["mult"] = sd.cache[:mult][sid]
            reacsys["info"]["chg"] = sd.cache[:charge][sid]
            reacsys["info"]["n_species"] = 1
        end
        @info "Assembled reactant system."
        if sum(rd.stoic_prods[i]) > 1
            prod_sids = []
            for (j, sid) in enumerate(rd.id_prods[i])
                for _ in 1:rd.stoic_prods[i][j]
                    push!(prod_sids, sid)
                end
            end
            prodsys = autode_NCI_conformer_search(sd, prod_sids; name="prodsys")
            if prodsys["info"]["chg"] != reacsys["info"]["chg"]
                throw(ErrorException("Charge not conserved in reaction $i: chg(R) = $(reacsys["info"]["chg"]), chg(P) = $(prodsys["info"]["chg"])"))
            end
            prodsys["info"]["n_species"] = length(prod_sids)
            prodsys_smi = join([sd.toStr[sid] for sid in prod_sids], ".")
            formal_charges = get_formal_charges(atom_map_smiles(prodsys, prodsys_smi))
            safe_geomopt!(prodsys, calc.calc_builder; calcdir=nebdir, mult=prodsys["info"]["mult"], 
                          chg=prodsys["info"]["chg"], formal_charges=formal_charges, maxiters=calc.maxiters)
        else
            sid = rd.id_prods[i][1]
            if sd.cache[:charge][sid] != reacsys["info"]["chg"]
                throw(ErrorException("Charge not conserved in reaction $i: chg(R) = $(reacsys["info"]["chg"]), chg(P) = $(sd.cache[:charge][sid])"))
            end
            prodsys = sd.xyz[sid]
            prodsys["info"]["mult"] = sd.cache[:mult][sid]
            prodsys["info"]["chg"] = sd.cache[:charge][sid]
            prodsys["info"]["n_species"] = 1
        end
        @info "Assembled product system."

        # Atom map endpoints.
        # Also obtain new formal charge arrays for remapped
        # endpoints. 
        reac_map, prod_map = split(rd.mapped_rxns[i], ">>")
        reac_map, prod_map = string(reac_map), string(prod_map)
        reacsys_mapped = atom_map_frame(reac_map, reacsys)
        reacsys_mapped["info"]["formal_charges"] = get_formal_charges(reac_map)
        prodsys_mapped = atom_map_frame(prod_map, prodsys)
        prodsys_mapped["info"]["formal_charges"] = get_formal_charges(prod_map)
        @info "Remapped atom indices."

        # Kabsch fit product system onto reactant system.
        kabsch_fit!(prodsys_mapped, reacsys_mapped)
        @info "Completed Kabsch fit of product system onto reactant system."

        # Interpolate and run NEB.
        images, conv = neb(reacsys_mapped, prodsys_mapped, calc; calcdir=nebdir)
        # Save to caches.
        # If unconverged and removal is requested, push blank
        # entries to caches so splice! still works at the end.
        if conv || !(calc.remove_unconverged)
            ts = highest_energy_frame(images)
            rxn_mult = get_rxn_mult(reacsys_mapped, prodsys_mapped)

            push!(calc.ts_cache[:xyz], ts)
            push!(calc.ts_cache[:mult], rxn_mult)
            push!(calc.ts_cache[:charge], prodsys_mapped["info"]["chg"])
            ts_sym, ts_geom = autode_frame_symmetry(ts; mult=rxn_mult, chg=prodsys_mapped["info"]["chg"])
            push!(calc.ts_cache[:symmetry], ts_sym)
            push!(calc.ts_cache[:geometry], ts_geom)
            rd.dH[i] = (prodsys_mapped["info"]["energy_ASE"] - reacsys_mapped["info"]["energy_ASE"]) * Constants.eV_to_kcal_per_mol
            push!(calc.ts_cache[:reacsys_energies], reacsys_mapped["info"]["energy_ASE"])
            push!(calc.ts_cache[:prodsys_energies], prodsys_mapped["info"]["energy_ASE"])

            # Run individual vibrational analyses on reactants, products
            # and TS.
            for rid in rd.id_reacs[i]
                calc_species_vibrations!(sd, rid, calc.calc_builder; calcdir=nebdir)
            end
            for pid in rd.id_prods[i]
                calc_species_vibrations!(sd, pid, calc.calc_builder; calcdir=nebdir)
            end
            calc_ts_vibrations!(calc.ts_cache, i, calc.calc_builder; calcdir=nebdir)
            @info "Completed vibrational analysis.\n"
        else
            push!(calc.ts_cache[:xyz], Dict{String, Any}())
            push!(calc.ts_cache[:vib_energies], [0.0+0.0im])
            push!(calc.ts_cache[:symmetry], -1)
            push!(calc.ts_cache[:geometry], -1)
            push!(calc.ts_cache[:reacsys_energies], 0.0)
            push!(calc.ts_cache[:prodsys_energies], 0.0)
            push!(calc.ts_cache[:mult], 0)
            push!(calc.ts_cache[:charge], 0)
        end

        calc.cached_rids[i] = true
        cd(currdir)
    end

    # Check all unconverged reactions for converged reverse pairs.
    # Copy the reverse TS if a converged pair is found.
    unconverged_rids = findall(==(-1), calc.ts_cache[:symmetry])
    for i in unconverged_rids
        reverse_rhash = get_reverse_rhash(sd, rd, i)
        if reverse_rhash in rd.rhash
            reverse_idx = findfirst(==(reverse_rhash), rd.rhash)
            if calc.ts_cache[:symmetry][reverse_idx] > -1 
                calc.ts_cache[:xyz][i] = calc.ts_cache[:xyz][reverse_idx]
                calc.ts_cache[:vib_energies][i] = calc.ts_cache[:vib_energies][reverse_idx]
                calc.ts_cache[:symmetry][i] = calc.ts_cache[:symmetry][reverse_idx]
                calc.ts_cache[:geometry][i] = calc.ts_cache[:geometry][reverse_idx]
                calc.ts_cache[:reacsys_energies][i] = calc.ts_cache[:prodsys_energies][reverse_idx]
                calc.ts_cache[:prodsys_energies][i] = calc.ts_cache[:reacsys_energies][reverse_idx]
                calc.ts_cache[:mult][i] = calc.ts_cache[:mult][reverse_idx]
                calc.ts_cache[:charge][i] = calc.ts_cache[:charge][reverse_idx]
                @info "Found a converged reverse reaction for RID $i, copied TS data."
            end
        end
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
        sd.cache[:mult][sid],
        sd.cache[:vib_energies][sid],
        T, P
    )
end

function get_entropy(ts_cache::Dict{Symbol, Any}, rid, mass, T, P)
    return get_entropy(
        mass,
        ts_cache[:xyz][rid]["info"]["inertias"],
        ts_cache[:geometry][rid],
        ts_cache[:symmetry][rid],
        ts_cache[:mult][rid],
        ts_cache[:vib_energies][rid], 
        T, P
    )
end

function get_entropy(mass, inertias, geometry, symmetry, mult, vib_energies, T, P)
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
    S += ASEConstants.kB * log10(mult)

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
function (calc::ASENEBCalculator{uType, tType})(; T::Number, P::Number) where {uType, tType}
    throw(ErrorException(""))
end

# Dispatched without k_max awareness.
function (calc::ASENEBCalculator{Nothing, tType})(; T::Number, P::Number) where {tType}
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

        S_ts = get_entropy(calc.ts_cache, rid, mass_ts, T, P)
        H_ts = get_enthalpy(calc.ts_cache, rid, T)

        dS[rid] = S_ts - S_reacs
        dH[rid] = H_ts - H_reacs
    end

    k = Constants.k_b*T/Constants.h .* exp.(dS/Constants.R) .* exp.(-dH/(Constants.R*T))
    return k
end

function Kinetica.has_conditions(::ASENEBCalculator, symbols::Vector{Symbol})
    return all([sym in [:T, :P] for sym in symbols])
end

Kinetica.allows_continuous(::ASENEBCalculator) = false