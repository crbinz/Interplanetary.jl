
export get_prop_mass,
       hohmann_circular_coplanar,
       hohmann_circular,
       hohmann_circle2circle,
       get_launch_asymptote,
       lambert,
       lambert_target,
       lambert_transfer

function get_prop_mass( m0::Float64,
                        dV::Float64,
                        Isp::Float64 ) 
    # dV in m/s, Isp in s, m0 in kg
    return m0 * (1. - exp(-dV/(9.81*Isp)))
end

"Hohmann transfer between two circular orbits"
function hohmann_circle2circle( a1::Float64,
                                a2::Float64,
                                mu::Float64 )
    aT = (a1 + a2)/2.
    vc1 = getv( mu, a1, a1 )
    vc2 = getv( mu, a2, a2 )
    
    vT1 = getv( mu, a1, aT )
    vT2 = getv( mu, a2, aT )

    dV1 = abs(vc1 - vT1)
    dV2 = abs(vc2 - vT2)

    return (aT, dV1, dV2)
end

function get_launch_asymptote( v_INF_HCI::Vector{Float64},
                               epoch::Float64 )
    T = get_hci2eci_matrix( epoch )
    v_INF_ECI = T*v_INF_HCI
    α = atan2(v_INF_ECI[2],v_INF_ECI[1])
    δ = asin(v_INF_ECI[3]/norm(v_INF_ECI))
    return (α,δ)
end

"""
Compute the relevant quantities for a Hohmann transfer between objects
orbiting the Sun, assuming circular, coplanar orbits also assuming
circular parking orbit
"""
function hohmann_circular_coplanar( muDepart::Float64,
                                    aDepart::Float64,
                                    rPark::Float64,
                                    muDest::Float64,
                                    aDest::Float64,
                                    aCapture::Float64,
                                    rCapture::Float64 )
    
    # v_INF_depart: v_infinity required wrt departure planet
    # v_INF_arrive: v_infinity wrt destination planet
    (aTrans, v_INF_depart, v_INF_arrive) = hohmann_circle2circle(aDepart, aDest, muSun)
    
    ### Departure

    # departure velocity wrt Sun:
    v_SC_S_depart = getv(muSun, aDepart, aTrans)
    
    # delta V required:
    dV_depart = getvhyp(muDepart, rPark, v_INF_depart) -
                sqrt(muDepart/rPark)
    
    ### Arrival
    
    # arrival velocity wrt Sun:
    v_SC_S_arrive = getv(muSun, aDest, aTrans)
    
    
    # delta V required for capture into orbit specified by aCapture and rCapture
    dV_cap = getArrivalDV( muDest, rCapture, aCapture, v_INF_arrive )
    
    TOF = getPeriod(muSun, aTrans)/2. # seconds
    
    return (v_INF_depart, dV_depart, v_INF_arrive, dV_cap, TOF)
end

function getArrivalDV( muDest::Float64,
                      rCapture::Float64,
                      aCapture::Float64,
                      v_INF_arrive::Float64 )
    return abs( getvhyp(muDest, rCapture, v_INF_arrive) -
               getv(muDest, rCapture, aCapture) )
end

"""
Compute the relevant quantities for a Hohmann transfer between objects
orbiting the Sun, assuming circular non-coplanar orbits
"""
function hohmann_circular( muDepart::Float64,
                           aDepart::Float64,
                           rPark::Float64,
                           muDest::Float64,
                           aDest::Float64,
                           dInc::Float64,
                           aCapture::Float64,
                           rCapture::Float64 )

    aTrans = (aDepart + aDest)/2
    
    # Departure

    # departure velocity wrt Sun:
    v_SC_S_depart = getv(muSun, aDepart, aTrans)
    v_depart_S = sqrt(muSun/aDepart)

    # arrival velocity wrt Sun:
    v_SC_S_arrive = getv(muSun, aDest, aTrans)
    v_dest_S = sqrt(muSun/aDest)

    # v_infinity required 
    v_INF_depart(alpha) = sqrt(v_SC_S_depart^2 + v_depart_S^2 - 
                        2*v_SC_S_depart*v_depart_S*cos(alpha))
    v_INF_arrive(alpha) = sqrt(v_SC_S_arrive^2 + v_dest_S^2 -
                               2*v_SC_S_arrive*v_dest_S*cos(dInc - alpha))

    # numerically solve for the optimal plane change split
    f(beta) = v_SC_S_depart*v_depart_S*sin(beta)/v_INF_depart(beta) - 
                v_dest_S*v_SC_S_arrive*sin(dInc - beta)/v_INF_arrive(beta)
    alpha1 = fzero(f,0.,dInc)

    # delta V required for departure:
    dV_depart = getvhyp(muDepart, rPark, v_INF_depart(alpha1))-
                sqrt(muDepart/rPark)
    
    # delta V required for capture into orbit specified by aCapture and rCapture
    dV_cap = getArrivalDV( muDest, rCapture, aCapture, v_INF_arrive(alpha1) )
    
    TOF = getPeriod(muSun, aTrans)/2. #seconds
    
    return (v_INF_depart(alpha1), dV_depart, v_INF_arrive(alpha1), dV_cap, alpha1, TOF)
end

"""
Translation of lambert_target.m

Original Matlab code provided by Brent Barbee
"""
function lambert_target(mu::Float64,
                        TOF::Float64,
                        r1::Vector{Float64},
                        r2::Vector{Float64},
                        v1::Vector{Float64},
                        v2::Vector{Float64})

    uu = cross(r1, r2)
    
    uh = (1/norm(uu))*uu
    
    (vi1, vf1) = lambert(mu, TOF, r1, r2, uh)

    tot_dv1 = norm(v1-vi1) + norm(vf1-v2)

    uu = cross(r2, r1)
    
    uh = (1/norm(uu))*uu
    
    (vi2, vf2) = lambert(mu, TOF, r1, r2, uh)

    tot_dv2 = norm(v1-vi2) + norm(vf2-v2)

    if tot_dv1 < tot_dv2

        vi = vi1
        vf = vf1

    else

        vi = vi2
        vf = vf2

    end
    return (vi, vf)
end

"""
Title     : Lambert                                                    
Purpose   : Given a conic gravity field, compute the required velocity 
            to transfer between initial and final inertial position    
            vectors in a specified transfer time interval.             
Inputs    : mu  - (km^3/s^2) Gravitational constant                    
            dt  - (s)        Transfer time interval                    
            ri  - (km)       Init. inertial position vector            
            rf  - (km)       Final inertial position vector            
            uh  - (na)       Unit angular momentum vector              
Outputs   : vi  - (km/s)     Init. inertial velocity vector (required) 
            vf  - (km/s)     Final inertial velocity vector            
Comments  : 1) This routine does not have a multi-revolution capability
            2) The unit angular momentum vector determines the plane,  
               and direction, of the transfer                          
Reference : Shepperd, S.W., unpublished notes.                         
Original author: Brent Barbee
"""
function lambert(mu::Float64,
                 dt::Float64,
                 ri::Vector{Float64},
                 rf::Vector{Float64},
                 uh::Vector{Float64})

    ### Initialization Iteration Loop 

    r0 = ri - dot(ri, uh) * uh
    r1 = rf - dot(rf, uh) * uh

    if abs(dot(uh, uh) - 1) > 10000 * eps()
        error("Lambert : Aborting : Normal vector not unitary\n")
    end

    if abs(dot(ri, uh)) > 0.00001
        error("Lambert : Warning : Init. position not in plane\n")
    end

    if abs(dot(rf, uh)) > 0.00001
        error("Lambert : Warning : Final position not in plane\n")
    end

    imax      = 20
    maxu      = typemax(Int64)
    tolerance = 1000 * eps()

    m0 = norm(r0)
    m1 = norm(r1)
    cc = norm(r1 - r0)
    ss = (m0 + m1 + cc) / 2

    rc = ss / 2
    vc = sqrt(mu / rc)
    wc = vc / rc

    k1 = sqrt(abs((ss - m0) * (ss - m1)))
    k  = dot(cross(r0, r1), uh)

    if  (k1 * k1) > (ss * (ss - cc))
        k  = k / (2 * k1 * ss)
        k2 = 1 - k * k
    else
        k2 = cc / ss
        if k >= 0
            k  = +sqrt(1 - k2)
        else
            k  = -sqrt(1 - k2)
        end
    end

    tdesired = wc * dt

    if tdesired > 4 * (1 - k * k * k) / 3
        u = 0 
        umin = -1 
        umax = 1
    else
        u = 1 
        umin = 1 
        umax = maxu
    end

    ### Iterate Until Transfer Time Matches 

    # initialize values outside of loop
    q = 0.
    y = 0.
    h1 = 0.
    h0 = 0.
    uu = 0.
    pp = 0.
    dt = 0.
    slope = 0.
    terror = 0.
    uold = 0.0
    dtold = 0.0

    for i=1:imax 

        q = k * u
        y = sqrt(q * q + k2)

        if q <= 0 
            h1 = y - q
        else
            h1 = k2 / (y + q)
        end

        h0 = k + u * h1
        q  = (1 - abs(h0)) / 2
        uu = (16 / 15) * h1 ^ 5 * cfgauss(5, q)

        if h0 < 0
            pp = 2 * pi / sqrt(1 - u * u)^5
            uu = pp - uu
        end

        dt     = 4 * h1 * (k + h0 * h1 * h1 / 3) + (1 - u * u) * uu
        slope  = 3 * u * uu - 4 * (h1 / y) * (k * k + (k * h0 + h1 * h1) * h1 * h1)
        terror = tdesired - dt

        if abs(terror) < wc
            if abs(terror) < tolerance * abs(tdesired )
                break
            end
            
            if abs(terror) < tolerance * abs(slope * u)
                break
            end
        end

        if (i > 1) && (u  ==  uold)
            break
        end

        if (i > 1) && (dt == dtold)
            break
        end

        uold  = u
        dtold = dt
        ustep = terror / slope

        if ustep > 0
            umin = u
            u    = u + ustep
            if u > umax
                u = (umin + umax) / 2
            end
        else
            umax = u
            u    = u + ustep
            if u < umin
                u = (umin + umax) / 2
            end
        end

    end

    ### Final Computations 

    h1 = h1 / vc
    uu = uu / vc^5

    h  = k1 / h1
    n0 = +(k * ss - m0 * h0) / h1
    n1 = -(k * ss - m1 * h0) / h1
    v0 = (n0 * r0 + h * cross(uh, r0)) / (m0 * m0)
    v1 = (n1 * r1 + h * cross(uh, r1)) / (m1 * m1)

    vi  = v0
    vf  = v1

    return (vi, vf)
end

"""
Title     : CF Gauss                                             
Purpose   : Compute the Gaussian continued fraction.             
Inputs    : i - (na) Order of continued fraction                 
            q - (na) Argument of continued fraction              
Outputs   : g - (na) Continued fraction result                   
Comments  : The continued fraction is defined as the ratio of two
            hypergeometric functions as follows:                 
                                                                 
                                   F(i,1,1+(i/2),q)              
                g = G(i,0,i/2,q) = ----------------              
                                     F(i,0,i/2,q)                
                                                                 
Exception : Illegal argument                                     
Reference : Shepperd, S.W., "Universal Keplerian State Transition
            Matrix," Celestial Mechanics, Vol. 35, 1985.         
Original author: Brent Barbee
"""
function cfgauss(i, q)

    if q > 0.5 + 16 * eps()
        error("CF Gauss : Illegal argument\n\n")
    end

    g = 1
    r = 1
    s = 1
    n = 0
    l = i - 2
    d = i * (i - 2)
    k = 1 - 2 * i

    while(true)

        k    = -k
        l    = l + 2
        d    = d + 4 * l
        n    = n + (1 + k) * l
        r    = d / (d - n * r * q)
        s    = (r - 1) * s
        gold = g
        g    = gold + s
        
        if g == gold
            break
        end

    end

    return g
end

    
"""
Wrapper to get delta V values for a Lambert-targeted transfer
assumes Sun as the central body

Units are km and s
"""
function lambert_transfer(muDepart::Float64, 
                         TOF::Float64,
                         launchEpoch::Float64,
                         rPark::Float64, 
                         x1_HCI::Vector{Float64}, x2_HCI::Vector{Float64},
                         muDest::Float64 = 0., 
                         rCapture::Float64 = 1000.,
                         aCapture::Float64 = 1000. )
    (v_SC_S_depart_HCI, v_SC_S_arrive_HCI) = lambert_target( muSun,
                                                             TOF,
                                                             x1_HCI[1:3], x2_HCI[1:3],
                                                             x1_HCI[4:6], x2_HCI[4:6] )
    # departure:
    v_INF_E_HCI = v_SC_S_depart_HCI - x1_HCI[4:6]
    C3_depart = norm(v_INF_E_HCI)^2
    dV_depart = abs(getvhyp(muEarth, rPark, norm(v_INF_E_HCI)) - getv(muEarth, rPark, rPark))
    
    (RLA,DLA) = get_launch_asymptote( v_INF_E_HCI, launchEpoch )
    
    # arrival:
    v_INF_arrive = norm(v_SC_S_arrive_HCI - x2_HCI[4:6])
    dV_arrive = abs(getArrivalDV( muDest, rCapture, aCapture, v_INF_arrive ))
    
    return (dV_depart, dV_arrive, C3_depart, RLA, DLA)
end
