# Module containing code for Interplanetary Navigation and Guidance class
# Christopher Binz
# University of Maryland

module Interplanetary

using Sundials, 
      Roots

export muSun, muEarth, AU, RE

const muSun   = 1.3271244004193938e11 # km^3/s^2
const muEarth = 3.98600440e5          # km^3/s^2
const AU      = 1.495978707e8         # km
const RE      = 6378.137              # km

include("utilities.jl")
include("frames.jl")
include("dynamics.jl")
include("transfers.jl")

end # end module

