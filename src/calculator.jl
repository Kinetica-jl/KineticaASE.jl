"""
"""
mutable struct ASENEBCalculator{kmType, uType, tType} <: Kinetica.AbstractKineticCalculator
    atoms_calc::Py
    neb_k
    climb::Bool
    climb_ftol
    maxiters::Int
    interpolation::Symbol
    k_max::kmType
    t_unit::String
    t_mult::tType
    cached_rids::Vector{Bool}
    cache::Dict{Symbol, Any}
end

"""
    ASENEBCalculator(atoms_calc[, neb_k, climb, climb_ftol, maxiters, interpolation, k_max, t_unit])

Outer constructor method for ASE-driven NEB-based kinetic calculator.
"""
function ASENEBCalculator(atoms_calc::Py; neb_k=0.1, climb::Bool=true, climb_ftol=0.04,
                          maxiters=200, interpolation::Symbol=:IDPP,
                          k_max::Union{Nothing, uType}=nothing, t_unit::String="s"
                         ) where {uType <: AbstractFloat}
    
    t_mult = tconvert(t_unit, "s")
    cached_rids = Vector{Bool}()
    cache = Dict{Symbol, Any}(
        :ts_atoms => Vector{Py}(),
        :reac_energies => Vector{Float64}(),
        :ts_energies => Vector{Float64}(),
        :geometries => Vector{String}(),
        :spins => Vector{String}(),
        :symmetries => Vector{Int}
    )
    return ASENEBCalculator(atoms_calc, neb_k, climb, climb_ftol, maxiters, 
                            interpolation, k_max, t_unit, t_mult, cached_rids,
                            cache)
end

function setup_network!(sd::SpeciesData, rd::RxData, calc::ASENEBCalculator)
    if length(calc.cached_rids) == 0
        calc.cached_rids = [false for _ in 1:rd.nr]
    elseif length(calc.cached_rids) < rd.nr
        calc.cached_rids = cat(calc.cached_rids, [false for _ in length(calc.cached_rids)+1:rd.nr]; dims=1)
    end

    for i in 1:rd.nr
        if calc.cached_rids[i] continue end
        reacsys = system_from_mols([sd.xyz[sid] for sid in rd.id_reacs[i]])
        prodsys = system_from_mols([sd.xyz[sid] for sid in rd.id_prods[i]])
        # Somehow create preconditioned systems and interpolate.
        # Run NEB.
        NEB = aseneb.NEB()
        # Run individual vibrational analyses on reac and ts fragments.
        # Calculate geometry, spin and symmetry of all fragments.
        #  - Consider storing these in their own SpeciesData?
        # Store everything in cache.
    end
end

function Base.splice!(calc::ASENEBCalculator, rids::Vector{Int})
    for key in keys(calc.cache)
        splice!(calc.cache[key], rids)
    end
    splice!(calc.cached_rids, rids)
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
    throw(ErrorException(""))
end

function has_conditions(::ASENEBCalculator, symbols::Vector{Symbol})
    return all([sym in [:T, :P] for sym in symbols])
end

allows_continuous(::ASENEBCalculator) = false