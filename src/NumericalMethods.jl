module NumericalMethods

include("derivative.jl")
include("ch5.jl")
include("ch6.jl")
include("ch9.jl")

export derivative
export bisect, false_position
export fixed_point, newton_raphson
export naive_gauss

end # module
