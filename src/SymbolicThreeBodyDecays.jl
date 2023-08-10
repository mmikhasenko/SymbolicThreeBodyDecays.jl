module SymbolicThreeBodyDecays
__precompile__(false)

using ThreeBodyDecay
import ThreeBodyDecay: cosθij, cosζ
import ThreeBodyDecay: MassTuple, MandestamTuple, WignerRotation, wr
# 
import ThreeBodyDecay.PartialWaveFunctions: wignerd_doublearg, CG_doublearg

using SymPy
import SymPy.PyCall

quantum_spin = PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
SymPy.import_from(quantum_spin, (:WignerD, :CG,), typ=:Any)

export StickySymTuple
export cosHold
include("dispatch.jl")

end # module SymbolicThreeBodyDecays
