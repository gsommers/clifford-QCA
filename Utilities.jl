#!/usr/licensed/julia/1.7/bin/julia

"""
General utility functions for data analysis and manipulation
"""
module Utilities

using JLD2
using LsqFit

export save_data, fit_linear, fit_function, get_binaries, get_base

lin(x,params)=(x.*params[1]).+ params[2]

lin_origin(x, params)=lin(x,[params,0])

"""
Fit ydata to the function `f` of xdata. 
Parameters
  - `f`: fitting function
  - `xdata`, `ydata`: x and y data points, respectively
  - `sigma`: uncertainties on y data points
  - `p0`: initial guess for parameters, i.e. y = f(xdata, p0)
  - `use_sigma::Bool` (kwarg): whether to include yerror in weighting the least square regression (default true)
  - `thresh` (kwarg): threshold for uncertainties: If `use_sigma` is `true`, but some uncertainties fall below `thresh`, reset them to `thresh` or they'll screw up the weights on the fits.
  - `print::Bool` (kwarg): whether to print warning about uncertainties below threshold.
Returns: best fit params, stderr, residuals
"""
function fit_function(f, xdata, ydata, sigma, p0; use_sigma::Bool=true, thresh=1e-6, print::Bool=false)
    if sigma==zeros(length(sigma)) || !use_sigma
        fitted = curve_fit(f, xdata, ydata, p0)
    else
	# these are bad for fitting
	zero_errs = findall(el->el<thresh,sigma)
	if !isnothing(zero_errs)
	    print && println("Warning: zero uncertainty at $(zero_errs)")
	    repl_sigma = copy(sigma)
	    repl_sigma[zero_errs] .= thresh
	else
	    repl_sigma = sigma
	end
	fitted = curve_fit(f,xdata, ydata, 1 ./repl_sigma.^2, p0)
    end
    # redefine resid to match what I would intuitively want
    resid = f(xdata, fitted.param) .- ydata
    fitted.param, stderror(fitted), resid
end

"""
Fit ydata to the linear function.
Parameters
  - `xdata`, `ydata`: x and y data points, respectively
  - `sigma`: uncertainties on y data points
  - `p0` (kwarg): initial guess for slope and intercept, i.e. y = p0[1] * xdata + p0[2]
  - `use_sigma::Bool` (kwarg): whether to include yerror in weighting the least square regression (default true)
  - `thresh` (kwarg): threshold for uncertainties: If `use_sigma` is `true`, but some uncertainties fall below `thresh`, reset them to `thresh` or they'll screw up the weights on the fits.
  - `print::Bool` (kwarg): whether to print warning about uncertainties below threshold.
Returns: best fit params (slope and intecept), stderr, residuals
"""
function fit_linear(xdata, ydata, sigma; p0 = [1.,0.],use_sigma::Bool=true,
	 thresh=1e-6, print::Bool=false)
    fit_function(lin, xdata, ydata, sigma, p0; use_sigma=use_sigma, thresh=thresh, print=print)
end

# returns: slope, stderr, residuals
function fit_proportional(xdata, ydata, sigma; p0 = [1.],use_sigma::Bool=true)
    fit_function(lin_origin, xdata, ydata, sigma, p0; use_sigma=use_sigma)
end

"""
Save data to .jld2 file with the name `name`. Optionally print message.
"""
function save_data(data, name;print::Bool=false)
    filename = "$(name).jld2"
    print && println("Saving to $(filename)")
    save(filename, data)
end

"""
Convert integer `num` into a boolean array with `places` number of bits.
"""
function get_binaries(num, places)
    Bool.(get_base(num, places, 2))
end

"""
Convert integer into an array of digits in a given base
""" 
function get_base(num, places, b)
    num < b^places ? nothing : throw(ArgumentError("Not enough places: $(num)>=b^$(places)"))
    bins = zeros(Int, places)
    idx = places;
    while num>0
        bins[idx] = num%b;
        num÷=b;
	idx-=1
    end
    bins
end

end