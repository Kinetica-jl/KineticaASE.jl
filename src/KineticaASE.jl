"""
KineticaASE.jl

UK Ministry of Defence © Crown Owned Copyright 2024/AWE
"""
module KineticaASE

using Kinetica
Constants = Kinetica.Constants

using Logging
using PythonCall
using ExtXYZ
using Glob
using StableHashTraits
using BSON

const version = VersionNumber(0, 1, 0)

const np = PythonCall.pynew()
const ase = PythonCall.pynew()
const aseopt = PythonCall.pynew()
const aseneb = PythonCall.pynew()
const aseio = PythonCall.pynew()
const asevib = PythonCall.pynew()
const asethermo = PythonCall.pynew()
const ade = PythonCall.pynew()
const rmsd = PythonCall.pynew()
function __init__()
    PythonCall.pycopy!(np, pyimport("numpy"))
    PythonCall.pycopy!(ase, pyimport("ase"))
    PythonCall.pycopy!(aseopt, pyimport("ase.optimize"))
    PythonCall.pycopy!(aseneb, pyimport("ase.neb"))
    PythonCall.pycopy!(aseio, pyimport("ase.io"))
    PythonCall.pycopy!(asevib, pyimport("ase.vibrations"))
    PythonCall.pycopy!(asethermo, pyimport("ase.thermochemistry"))
    PythonCall.pycopy!(ade, pyimport("autode"))
    PythonCall.pycopy!(rmsd, pyimport("rmsd"))
end

include("constants.jl")
using .ASEConstants

include("conversion.jl")
export frame_to_atoms, atoms_to_frame

include("calculator.jl")
export ASENEBCalculator
export calculate_entropy_enthalpy

include("neb.jl")

include("optimise.jl")

include("vibrations.jl")

include("autode.jl")

include("builders.jl")
export EMTBuilder, NWChemDFTBuilder

include("io.jl")

end