# dynamics functions, including the restricted three-body problem

export r3bp_hci_eom, twobody_eom, propagate_twobody, propagate_r3bp_hci

"""
Derivatives function for the restricted three-body problem in the
Heliocentric Inertial (HCI) reference frame xHCI is the stacked state
vector for the secondary planet + satellite: [xPlanet, xSat]
"""
function r3bp_hci_eom( t::Float64,
                       xHCI::Vector{Float64},
                       xdot,
                       mu2::Float64 = muEarth )
    
    derivatives = zeros(12)
    derivatives[1:3] = xHCI[4:6]   # Planet's velocity
    derivatives[7:9] = xHCI[10:12] # Satellite velocity
    
    rPlanet = xHCI[1:3]
    rSat    = xHCI[7:9]
    derivatives[4:6]   = -muSun/(norm(rPlanet)^3) .* rPlanet
    derivatives[10:12] = -muSun/(norm(rSat)^3) .* rSat -
                            muEarth/(norm(rSat - rPlanet)^3) .* (rSat - rPlanet)
    xdot[:]=derivatives
end

"Derivatives function for the two-body problem"
function twobody_eom( t::Float64,
                      x::Vector{Float64},
                      xdot,
                      mu::Float64 = muEarth )
    
    xdot[1:3] = x[4:6]
    xdot[4:6] = -mu/(norm(x[1:3])^3) * x[1:3]
end

"""
Propagate the R3BP in the HCI frame vectors are expected in the HCI
frame 

t should be an array of time values
"""
function propagate_r3bp_hci( x0Planet::Vector{Float64},
                             x0Sat::Vector{Float64},
                             t )

    return Sundials.cvode(r3bp_hci_eom,
                   vcat(x0Planet, x0Sat),
                   t, 
                   abstol=1e-12,reltol=1e-12)'
end

"Propagate the two-body problem"
function propagate_twobody( x0::Vector{Float64},
                            t )
    
    return Sundials.cvode(twobody_eom,
                          x0,
                          t,
                          abstol=1e-12,reltol=1e-12)'
end
