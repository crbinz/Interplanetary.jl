# This is an example which uses much of the functionality provided by
# Interplanetary.jl. It is based on a homework assignment, but that is
# no guarentee that it's correct.
#
# The objective is to explore optimal trajectories from Earth to
# asteroid 2015PDC - a fictional asteroid created for a hypothetical
# asteroid impact scenario held at the 2015 Planetary Defense
# Conference (see http://neo.jpl.nasa.gov/pdc15/). This is
# accomplished via a Pork Chop Plot - a contour plot showing the value
# of some objective function (mass at arrival, delta V, etc.) over a
# grid of Earth departure dates and time of flight values.
#
# This particular example searches for the minimum-delta V solution. I
# also included some of my experiments with using Julia's built-in
# parallel computing capabilities, because it seemed like a good place
# for parallelization.

function meshgrid( x::Vector{Float64}, y::Vector{Float64} )
    return ( repmat(x, 1, length(y)), 
             repmat(y',length(x),1) )
end

"""
This function creates a Pork Chop Plot given a range of departure
dates, time of flight for a given grid point (trajectory), and the
figure of merit value for a given trajectory.
"""
function createPCC( departure_date_range::Vector{Float64},
                   TOF_range::Vector{Float64},
                   FOM_array::Matrix{Float64}; levels = [],
                   FOM_label::ASCIIString = "")
    (depart_grid, TOF_grid) = meshgrid(departure_date_range, TOF_range)

    if isempty(levels)
        cplot = contourf(depart_grid, TOF_grid, FOM_array)
    else
        cplot = contourf(depart_grid, TOF_grid, FOM_array, levels)
    end
    
    # formatting
    ylocs,ylabels = yticks()
    xlocs,xlabels = xticks()
    cb = colorbar()
    cb[:set_label](FOM_label)

    xticks(xlocs,Date(julian2datetime(mjd2jd(xlocs))),rotation=30)
    #yticks(ylocs,Date(julian2datetime(mjd2jd(ylocs))))
    xlabel("Departure date")
    ylabel("TOF [days]")
    # need to set axis limits
end

# read ephemeris files for Earth and Asteroid 2015PDC
Earth_ephem = readdlm("data/Earth_2015-04-01_to_2022-09-03_step_1d.orb")
ast_2015_PDC_ephem = readdlm("data/ast_2015_PDC_2015-04-01_to_2022-09-03_step_1d.orb")

@everywhere begin # need to add cores before running this: addprocs(n)
using Dates, Interplanetary, PyPlot
departure_lb_mjd = jd2mjd( datetime2julian( DateTime("2015-10-01") ) )
departure_ub_mjd = jd2mjd( datetime2julian( DateTime("2022-02-26") ) )

# set up ranges for grid search
departure_step = 5. #days
TOF_step = 5. #days

departure_date_range = [departure_lb_mjd:departure_step:departure_ub_mjd] 
TOF_range = [90.:TOF_step:1500.]

# now loop over all combinations

#dV_total = Array(Float64,N,N); # for non-parallel
dV_total = SharedArray(Float64,length(departure_date_range),length(TOF_range))
                 
rPark = RE + 185. #parking orbit radius
end

@time @sync @parallel for i = 1:length(departure_date_range)
    for j = 1:length(TOF_range)
        dep_date = departure_date_range[i]
        arr_date = dep_date + TOF_range[j]
        if arr_date <= ast_2015_PDC_ephem[end,1]
            (dv1,dv2,c3,rla,dla) = lambert_transfer(muEarth,
                                        TOF_range[j]*86400.,
                                        mjd2jd(departure_date_range[i]),
                                        rPark,
                                        vec(Earth_ephem[Earth_ephem[:,1].==dep_date,2:7]),
                                        vec(ast_2015_PDC_ephem[ast_2015_PDC_ephem[:,1].==arr_date,2:7]))
            dV_total[i,j] = dv1+dv2
        end#if
    end#for
end#for

# create the contour plot
# sdata() makes a SharedArray into a regular array
createPCC(departure_date_range,
          TOF_range,sdata(dV_total),
          levels = [5.:2:30],
          FOM_label = "Delta V [km/s]") 

# find the minimum dV solution
# linear index of minimum, excluding zeros:
ind = find(dV_total .== minimum(dV_total[dV_total.>0])) 

subs = ind2sub(size(dV_total),ind) # convert to array subscripts
opt_arrival_date = departure_date_range[subs[1][1]] + TOF_range[subs[2][1]]

# plot the optimal solution on the pork chop plot
plot(departure_date_range[subs[1][1]],TOF_range[subs[2][1]],"r*",markersize=10)

# rerun the Lambert targeter for the optimal solution
(dv1,dv2,c3,rla,dla) = lambert_transfer(muEarth,
                            TOF_range[subs[2][1]]*86400.,
                            mjd2jd(departure_date_range[subs[1][1]]),
                            rPark,
                            vec(Earth_ephem[Earth_ephem[:,1].==departure_date_range[subs[1][1]],2:7]),
                            vec(ast_2015_PDC_ephem[ast_2015_PDC_ephem[:,1].==opt_arrival_date,2:7]))
                                       
                                       
# print results to screen
println("Optimal solution:")
println("Departure date: ", Date(julian2datetime(mjd2jd(departure_date_range[subs[1][1]]))))
println("TOF (days): ",TOF_range[subs[2][1]])
println("Earth departure ΔV [km/s]: ", dv1)
println("Earth departure C3 [km^2/s^2]: ", c3)
println("DLA (deg): ", dla*180./pi)
println("Arrival ΔV [km/s]: ", dv2)
println("Total solution ΔV [km/s]: ", dV_total[subs...][1])
