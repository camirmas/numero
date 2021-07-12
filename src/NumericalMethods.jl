module NumericalMethods

include("derivative.jl")
include("ch5.jl")
include("ch6.jl")

export derivative
export bisect, false_position
export fixed_point, newton_raphson

end # module
