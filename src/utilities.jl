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

function rv2kep( x_rv::Vector{Float64}, mu::Float64 = muEarth )
    # ref: GTDS Math Spec; Vallado, 2nd ed., p. 121
    # *NOTE*: undefined quantities return 0.0
    r = norm(x_rv[1:3])
    v = norm(x_rv[4:6])
    h = cross(x_rv[1:3],x_rv[4:6]) # angular momentum vector
    hmag = norm(h)                 # magnitude of angular momentum
    n = cross([0;0.;1],h)         
    
    energy = v^2/2. - mu/r
    
    #eccentricity vector
    evec = ((v^2 - mu/r)*x_rv[1:3] - dot(x_rv[1:3],x_rv[4:6])*x_rv[4:6])/mu
    ecc = norm(evec) #eccentricity
    
    # semimajor axis
    if abs(ecc - 1.) > eps()
        sma = -mu/(2.*energy)
    else
        sma = Inf
    end

    # inclination
    inc = acos(round(h[3]/hmag,10)) #round to avoid floating point issues
    
    # handle special cases
    if norm(n) < eps() # equatorial
        RAAN = 0.
    else
        RAAN = acos(n[1]/norm(n))
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
            AOP = acos(round(dot(n,evec)/(norm(n)*ecc),10)) # workaround of floating point issues
            if evec[3]<0.
                AOP = 2pi - AOP
            end
        end
        nu = acos(round(dot(evec,x_rv[1:3])/(ecc*r),10))
        if dot(x_rv[1:3],x_rv[4:6]) < 0.
            nu = 2pi - nu
        end
    end
    return [sma;ecc;inc;RAAN;AOP;nu]
end

"Wrapper for rv2kep to handle a vector of Cartesian vectors"
function rv2kep( x_rv::Matrix{Float64} )
   xkep = zeros(size(x_rv))
   for i = 1:size(x_rv,2)
       xkep[:,i] = rv2kep(x_rv[:,i])
   end
   return xkep
end

function get_hyperbolic_anomaly( r::Float64, a_hyp::Float64, e_hyp::Float64 )
    return acosh((1. - r/a_hyp)/e_hyp)
end

"""
Return the TOF in seconds given the SMA and eccentricity of a
hyperbola and r0 and r1
"""
function get_hyperbolic_tof( r0::Float64, r1::Float64, a_hyp::Float64, e_hyp::Float64 )
    H0 = get_hyperbolic_anomaly(r0,a_hyp,e_hyp)
    H1 = get_hyperbolic_anomaly(r1,a_hyp,e_hyp)
    
    return sqrt(-a_hyp^3/muEarth) * 
            (e_hyp * sinh(H1) - H1 - (e_hyp * sinh(H0) - H0))
end

"""
Calculate the Sphere of Influence radius for a body, given it's
properties and the properties of the larger body
"""
function get_soi( m_smaller, M_larger, R_CB )
    return R_CB*(m_smaller/M_larger)^(2/5.)
end

"""
Calculate the Sphere of Activity radius for a body, given it's
properties and the properties of the larger body
"""
function get_soa( m_smaller, M_larger, R_CB )
    return R_CB*(m_smaller/(3.*M_larger))^(1/3)
end

"Calculate elliptical orbit velocity"
function getv( mu::Float64, r::Float64, a::Float64 )
    return sqrt(mu*((2/r) - 1/a))
end

"Calculate hyperbolic orbit velocity"
function getvhyp( mu::Float64, r::Float64, vInf::Float64 )
    return sqrt(vInf^2 + 2*mu/r)
end

function getperiod( mu::Float64, a::Float64 )
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

function mjd2jd(mjd::Vector{Float64})
    jd = zeros(size(mjd))
    for i = 1:length(mjd)
        jd[i] = mjd2jd(mjd[i])
    end
    return jd
end

function jd2mjd(jd::Vector{Float64})
    mjd = zeros(size(jd))
    for i = 1:length(jd)
        mjd[i] = jd2mjd(jd[i])
    end
    return mjd 
end

"Calculate the approach phase angle (Sun - x_target - x)"
function get_approach_phase_angle( x::Vector{Float64}, x_target::Vector{Float64} )
    v_rel = (x[4:6] - x_target[4:6])./norm(x[4:6] - x_target[4:6])
    return acos( dot(v_rel, x_target[1:3]./norm(x_target[1:3]) ) )
end

"Calculate the Sun-Spacecraft-Earth angle"
function get_sse( r_SC_S::Vector{Float64}, r_E_S::Vector{Float64} )
    r_S_SC = -r_SC_S
    r_E_SC = r_E_S - r_SC_S
    return acos( dot(r_S_SC,r_E_SC)/(norm(r_S_SC)*norm(r_E_SC)) )
end

"Parse a HORIZONS ephemeris file"
function get_HORIZONS_ephem( f::IOStream )
    ephem = Array(Float64,0,7)
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
                ephem = vcat( ephem, [jdate_current eph_current] )
                x= readline(f)
            end
        end
    end
    return ephem
end  
