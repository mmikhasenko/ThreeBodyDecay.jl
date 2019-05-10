module ThreeBodyDecay

export σ3of1, σ1of2, σ2of3,
       cosθ12, cosθ23, cosθ31,
       cosθhat12,
       cosβ,
       λ,
       Kibble, Kibble23, Kibble12, Kibbl31

export ThreeBodySystem
export pthsq, pth

include("angle_functions.jl")
include("tbs_struct.jl")

end  # module ThreeBodyDecay
