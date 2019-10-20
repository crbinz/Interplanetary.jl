export rv2kep,
        get_hyperbolic_tof,
        get_hyperbolic_anomaly,
        get_soa,
        get_soi,
        getv,
        getvhyp,
        mjd2jd,
        jd2mjd,
        get_approach_phase_angle,
        get_HORIZONS_ephem,
        get_sse

"""
Transform a Cartesian ECI state vector to a set of
classical/Keplerian elements.

x_kep - [a,e,i,RAAN,AOP,nu], where:
    a    = semimajor axis [km]
    e    = eccentricity
    i    = inclination [rad]
    RAAN = right ascension of the ascending node [rad]
    AOP  = argument of perigee [rad]
    nu   = true anomaly [rad]

Reference: GTDS math spec

**Undefined quantities return as 0.0**
"""
function rv2kep(x_rv::AbstractVector{T}, mu = muEarth) where {T<:AbstractFloat}
    r = norm(x_rv[1:3])
    v = norm(x_rv[4:6])
    h = cross(x_rv[1:3],x_rv[4:6])
    n = cross([0;0.;1],h)

    sma = mu * r/(2*mu - r*v^2)
    evec = eccvec(x_rv)
    ecc = norm(evec)
    inc = angle_between(h, [0.0; 0.0; 1.0])

    if norm(n) < eps() # equatorial
        RAAN = 0.
    else
        RAAN = angle_between([1.0; 0.0; 0.0], n)
        if n[2] < 0.
            RAAN = 2pi - RAAN
        end
    end

    # circular
    if ecc < eps()
        AOP = 0.
        nu = 0.
    else
        if norm(n) < eps() # equatorial
            AOP = 0.0
        else
            AOP = angle_between(n, evec)
            if evec[3]<0.
                AOP = 2pi - AOP
            end
        end
        nu = angle_between(evec, x_rv[1:3])
        if dot(x_rv[1:3],x_rv[4:6]) < 0.
            nu = 2pi - nu
        end
    end
    return [sma;ecc;inc;RAAN;AOP;nu]
end

function get_hyperbolic_anomaly(r, a_hyp, e_hyp)
    return acosh((1. - r/a_hyp)/e_hyp)
end

"""
Return the TOF in seconds given the SMA and eccentricity of a
hyperbola and r0 and r1
"""
function get_hyperbolic_tof(r0, r1, a_hyp, e_hyp)
    H0 = get_hyperbolic_anomaly(r0,a_hyp,e_hyp)
    H1 = get_hyperbolic_anomaly(r1,a_hyp,e_hyp)
    
    return sqrt(-a_hyp^3/muEarth) * 
            (e_hyp * sinh(H1) - H1 - (e_hyp * sinh(H0) - H0))
end

"""
Calculate the Sphere of Influence radius for a body, given it's
properties and the properties of the larger body
"""
function get_soi(m_smaller, M_larger, R_CB)
    return R_CB*(m_smaller/M_larger)^(2/5.)
end

"""
Calculate the Sphere of Activity radius for a body, given it's
properties and the properties of the larger body
"""
function get_soa( m_smaller, M_larger, R_CB )
    return R_CB*(m_smaller/(3.0 * M_larger))^(1/3)
end

"Calculate elliptical orbit velocity"
function getv(mu, r, a)
    return sqrt(mu*((2/r) - 1/a))
end

"Calculate hyperbolic orbit velocity"
function getvhyp(mu, r, vInf)
    return sqrt(vInf^2 + 2*mu/r)
end

function getperiod(mu, a)
    return 2pi*sqrt(a^3/mu)
end

"Convert Julian date to Modified Julian date"
function jd2mjd(jd::Float64)
    return jd - 2400000.5
end

"Convert Modified Julian date to Julian date"
function mjd2jd(mjd::Float64)
    return mjd + 2400000.5
end

"Calculate the approach phase angle (Sun - x_target - x)"
function get_approach_phase_angle(x::Vector{T}, x_target::Vector{T}) where {T<:AbstractFloat}
    v_rel = (x[4:6] - x_target[4:6])./norm(x[4:6] - x_target[4:6])
    return acos( dot(v_rel, x_target[1:3]./norm(x_target[1:3]) ) )
end

"Calculate the Sun-Spacecraft-Earth angle"
function get_sse(r_SC_S::Vector{T}, r_E_S::Vector{T}) where {T<:AbstractFloat}
    r_S_SC = -r_SC_S
    r_E_SC = r_E_S - r_SC_S
    return acos( dot(r_S_SC,r_E_SC)/(norm(r_S_SC)*norm(r_E_SC)) )
end

"Parse a HORIZONS ephemeris file"
function get_HORIZONS_ephem(f::IOStream)
    ephem = Vector{Vector{Float64}}(undef,0)
    regex_jdate = r"\d{7}.\d{9}"
    #regex_datetime = r"\d{4}-[a-zA-Z]{3}-\d{2}" # not used
    #regex_time = r"\d{2}:\d{2}:\d{2}.\d+"       # not used
    regex_state = r"[-+]?\d+?.\d+[eE][+-]\d+.*$" #matches from the first match to the end of the line
    while !eof(f)
        x = readline(f)
        if contains(x,"\$\$SOE") # start of ephem
            x = readline(f)
            while !contains(x,"\$\$EOE") #end of ephemeris
                eph_current = float64(split(match(regex_state,x).match,',')[1:6])' #extra comma is off the end
                jdate_current = float64(match(regex_jdate,x).match)
                push!(ephem, vcat(jdate_current,eph_current))
                x= readline(f)
            end
        end
    end
    return hcat(ephem...)'
end  
