"""
KineticaASE.jl

UK Ministry of Defence © Crown Owned Copyright 2024/AWE
"""
module KineticaASE

using Kinetica
Constants = Kinetica.Constants

using Logging
using PythonCall

const version = VersionNumber(0, 1, 0)

const ase = PythonCall.pynew()
const aseneb = PythonCall.pynew()
const aseio = PythonCall.pynew()
const asevib = PythonCall.pynew()
const asethermo = PythonCall.pynew()
function __init__()
    PythonCall.pycopy!(ase, pyimport("ase"))
    PythonCall.pycopy!(aseneb, pyimport("ase.neb"))
    PythonCall.pycopy!(aseio, pyimport("ase.io"))
    PythonCall.pycopy!(asevib, pyimport("ase.vibrations"))
    PythonCall.pycopy!(asethermo, pyimport("ase.thermochemistry"))
end

end