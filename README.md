# KineticaASE.jl

A modular kinetic calculator addon for [Kinetica.jl](https://github.com/Kinetica-jl/Kinetica.jl), enabling calculation of temperature and pressure-dependent rates of reaction through the Atomic Simulation Environment [(ASE)](https://wiki.fysik.dtu.dk/ase/index.html).

This package makes use of ASE's implementation of the nudged elastic band (NEB) method to determine transition states of every reaction in a chemical reaction network (CRN), then characterises reactions through vibrational analysis to calculate their Gibbs free energies of activation $\Delta G^{\ddagger}$, which can be used to calculate reaction rates.