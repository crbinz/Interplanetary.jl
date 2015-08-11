# Reference frame transformations

export get_hci2eci_matrix

"""
Compute the rotation matrix for the HCI -> ECI transformation 

jd is the Julian date
"""
function get_hci2eci_matrix( jd::Float64 )
    
    T_JC = (jd - 2451545.0)/36525.

    eps = 0.409092802283074 - 0.000226966106587847 * T_JC - 
    2.8623399732707e-9 * T_JC^2 + 
    8.79645943005142e-9 * T_JC^3

    # Rotation matrix from HCI to ECI
    return [1. 0 0;
            0 cos(eps) -sin(eps);
            0 sin(eps) cos(eps)]
end
