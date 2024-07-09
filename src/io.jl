"""
"""
function save_asecalc(calc::ASENEBCalculator, saveto::String)
    savedict = Dict(
        :KineticaASEVersion => pkgversion(KineticaASE),
        :calcdir_head => calc.calcdir_head,
        :neb_k => calc.neb_k,
        :ftol => calc.ftol,
        :climb => calc.climb,
        :climb_ftol => calc.climb_ftol,
        :maxiters => calc.maxiters,
        :interpolation => calc.interpolation,
        :n_images => calc.n_images,
        :parallel => calc.parallel,
        :neb_optimiser => calc.neb_optimiser,
        :remove_unconverged => calc.remove_unconverged,
        :vibration_displacement => calc.vibration_displacement,
        :imaginary_freq_tol => calc.imaginary_freq_tol,
        :k_max => calc.k_max,
        :t_unit => calc.t_unit,
        :t_mult => calc.t_mult,
        :cached_rhashes => calc.cached_rhashes,

        :ts_cache => calc.ts_cache,
        :sd => Dict(
            :toInt => calc.sd.toInt,
            :n => calc.sd.n,
            :xyz => calc.sd.xyz,
            :cache => calc.sd.cache
        ),
        :rd => Dict(
            :nr => calc.rd.nr,
            :mapped_rxns => calc.rd.mapped_rxns,
            :id_reacs => calc.rd.id_reacs,
            :id_prods => calc.rd.id_prods,
            :stoic_reacs => calc.rd.stoic_reacs,
            :stoic_prods => calc.rd.stoic_prods,
            :dH => calc.rd.dH,
            :rhash => calc.rd.rhash
        )
    )

    bson(saveto, savedict)
end


"""
"""
function load_asecalc(calcfile::String, calc_builder)
    savedict = BSON.load(calcfile)
    if savedict[:KineticaASEVersion] < pkgversion(KineticaASE)
        @warn "Loaded calculator was made by a previous version of KineticaASE, check reconstruction is correct."
    elseif savedict[:KineticaASEVersion] > pkgversion(KineticaASE)
        @warn "Loaded calculator was made by a newer version of KineticaASE, check reconstruction is correct." 
    end

    sd_iType = typeof(savedict[:sd][:n])
    sd_toStr = Dict{sd_iType, String}(value => key for (key, value) in savedict[:sd][:toInt])
    sd = SpeciesData(
        savedict[:sd][:toInt],
        sd_toStr,
        savedict[:sd][:n],
        savedict[:sd][:xyz],
        savedict[:sd][:cache]
    )

    rd_iType = typeof(savedict[:rd][:nr])
    rd_mapped_rxns = get(savedict[:rd], :mapped_rxns, String[])
    if length(rd_mapped_rxns) == 0
        @warn "No reaction atom maps found in output."
    end
    rd = RxData(
        savedict[:rd][:nr], String[rxn for rxn in rd_mapped_rxns],
        Vector{rd_iType}[reac for reac in savedict[:rd][:id_reacs]], 
        Vector{rd_iType}[prod for prod in savedict[:rd][:id_prods]],
        Vector{rd_iType}[sreac for sreac in savedict[:rd][:stoic_reacs]], 
        Vector{rd_iType}[sprod for sprod in savedict[:rd][:stoic_prods]],
        savedict[:rd][:dH], 
        Vector{UInt8}[hash for hash in savedict[:rd][:rhash]]
    )

    # Work around Vector{Vector} type instability.
    rhashes = Vector{UInt8}[rh for rh in savedict[:cached_rhashes]]

    calc = ASENEBCalculator(
        calc_builder, 
        savedict[:calcdir_head],
        savedict[:neb_k],
        savedict[:ftol],
        savedict[:climb],
        savedict[:climb_ftol],
        savedict[:maxiters],
        savedict[:interpolation],
        savedict[:n_images],
        savedict[:parallel],
        savedict[:neb_optimiser],
        savedict[:remove_unconverged],
        savedict[:vibration_displacement],
        savedict[:imaginary_freq_tol],
        savedict[:k_max],
        savedict[:t_unit],
        savedict[:t_mult],
        rhashes,
        savedict[:ts_cache],
        sd,
        rd
    )
    return calc
end


"""
"""
function load_asecalc_data!(calc::ASENEBCalculator, calcfile::String)
    savedict = BSON.load(calcfile)
    if savedict[:KineticaASEVersion] < pkgversion(KineticaASE)
        @warn "Loaded calculator was made by a previous version of KineticaASE, check reconstruction is correct."
    elseif savedict[:KineticaASEVersion] > pkgversion(KineticaASE)
        @warn "Loaded calculator was made by a newer version of KineticaASE, check reconstruction is correct." 
    end

    sd_iType = typeof(savedict[:sd][:n])
    sd_toStr = Dict{sd_iType, String}(value => key for (key, value) in savedict[:sd][:toInt])
    sd = SpeciesData(
        savedict[:sd][:toInt],
        sd_toStr,
        savedict[:sd][:n],
        savedict[:sd][:xyz],
        savedict[:sd][:cache]
    )

    rd_iType = typeof(savedict[:rd][:nr])
    rd_mapped_rxns = get(savedict[:rd], :mapped_rxns, String[])
    if length(rd_mapped_rxns) == 0
        @warn "No reaction atom maps found in output."
    end
    rd = RxData(
        savedict[:rd][:nr], String[rxn for rxn in rd_mapped_rxns],
        Vector{rd_iType}[reac for reac in savedict[:rd][:id_reacs]], 
        Vector{rd_iType}[prod for prod in savedict[:rd][:id_prods]],
        Vector{rd_iType}[sreac for sreac in savedict[:rd][:stoic_reacs]], 
        Vector{rd_iType}[sprod for sprod in savedict[:rd][:stoic_prods]],
        savedict[:rd][:dH], 
        Vector{UInt8}[hash for hash in savedict[:rd][:rhash]]
    )

    calc.ts_cache = savedict[:ts_cache]
    calc.sd = sd
    calc.rd = rd
    return
end


"""
"""
function verify_sd(calc_sd::SpeciesData, real_sd::SpeciesData)
    if calc_sd.n == 0 return end
    if calc_sd.n > real_sd.n
        throw(ErrorException("Calculator SpeciesData is not a subset of input SpeciesData."))
    end

    mismatch_sids = []
    for i in 1:calc_sd.n
        if calc_sd.toStr[i] != real_sd.toStr[i]
            push!(mismatch_sids, i)
        end
    end
    if length(mismatch_sids) > 0
        throw(ErrorException("Calculator SpeciesData deviates from input SpeciesData at SIDs $(mismatch_sids)"))
    end
    return
end


"""
"""
function verify_rd(calc_rd::RxData, real_rd::RxData)
    if calc_rd.nr == 0 return end
    if calc_rd.nr > real_rd.nr
        throw(ErrorException("Calculator RxData is not a subset of input RxData."))
    end

    mismatch_rids = []
    for i in 1:calc_rd.nr
        if calc_rd.rhash[i] != real_rd.rhash[i]
            push!(mismatch_rids, i)
        end
    end
    if length(mismatch_rids) > 0
        throw(ErrorException("Calculator RxData deviates from input RxData at RIDs $(mismatch_rids)"))
    end
    return
end


"""
"""
function save_optgeom(frame::Dict{String, Any}, sym, geom, saveto::String)
    savedict = Dict(
        :frame => frame,
        :sym => sym,
        :geom => geom
    )
    bson(saveto, savedict)
end


"""
"""
function load_optgeom(savefile::String)
    savedict = BSON.load(savefile)
    savedict[:frame]["arrays"]["species"] = String[s for s in savedict[:frame]["arrays"]["species"]]
    return savedict
end



"""
"""
function save_endpoints(reacsys::Dict{String, Any}, prodsys::Dict{String, Any}, saveto::String)
    savedict = Dict(
        :reacsys => reacsys,
        :prodsys => prodsys
    )
    bson(saveto, savedict)
end


"""
"""
function load_endpoints(savefile::String)
    savedict = BSON.load(savefile)
    return savedict[:reacsys], savedict[:prodsys]
end


"""
"""
function save_tsdata(ts::Dict{String, Any}, conv, mult, sym, geom, chg, saveto::String)
    savedict = Dict(
        :conv => conv,
        :xyz => ts,
        :mult => mult,
        :charge => chg,
        :sym => sym,
        :geom => geom
    )
    bson(saveto, savedict)
end


"""
"""
function load_tsdata(savefile::String)
    savedict = BSON.load(savefile)
    return savedict
end


"""
"""
function save_vibdata(sd::SpeciesData, sids, ts_cache::Dict{Symbol, Any}, rid, saveto::String)
    savedict = Dict(
        :smis => [sd.toStr[sid] for sid in sids],
        :by_smi => Dict(sd.toStr[sid] => sd.cache[:vib_energies][sid] for sid in sids),
        :ts => ts_cache[:vib_energies][rid]
    )
    bson(saveto, savedict)
end


"""
"""
function load_vibdata(savefile::String)
    savedict = BSON.load(savefile)
    return savedict
end