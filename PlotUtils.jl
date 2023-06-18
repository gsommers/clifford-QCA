#!/usr/licensed/julia/1.7/bin/julia

"""
Utility functions for making and labeling plots
"""
module PlotUtils
using Statistics
using Plots
using Plots.PlotMeasures
using Utilities
using LaTeXStrings
using Formatting

export plot_suptitle, plot_layout, plot_fit_linear, get_legend_key

function get_legend_key(L; lab = "L")
    if !(typeof(L)<:Number)
        return L
    end
    if L%1==0
        L = Int(L)
    end
    if typeof(L)<:Int
       return LaTeXString(format(L"{}={:d}",lab,L));
    else
       return LaTeXString(format(L"{}=", lab)) * "$(round(L,digits=3))"
    end
end

"""
Add a title to the plot plt, which possibly consists of many panels. If a nonzero value of num_plots is passed in, automatically determine the plot size to use, based on how many subplots are in plt. Returns the plot p
"""
function plot_suptitle(plt, title; size=(900,750),num_plots::Int=0, fontsize=12, rel_size=(1,0.9))
     y = ones(3)
     plt_title = Plots.scatter(y, marker=0,markeralpha=0, 
     	annotations=(2, 1.5, Plots.text(title, fontsize)),axis=false, grid=false,legend=false, top_margin=50px)
     size = get_plot_size(use_size=size, num_plots=num_plots, rel_size=rel_size)
     p = Plots.plot(plt_title, plt, layout=Plots.grid(2,1,heights=[0.1*600/size[2],1-0.1*600/size[2]]), size = size)
 end

function plot_layout(all_plots, ncols; force_grid::Bool=true)
    if !force_grid || ncols==1 || ncols==length(all_plots) # don't make a grid
        p = Plots.plot(all_plots...)
	use_size = get_plot_size(num_plots=length(all_plots))
    else
	nrows = length(all_plots)÷ncols
        p = Plots.plot(all_plots...,layout=(nrows,ncols))
	use_size = (300*ncols,300*nrows)
    end
    return p, use_size
end

# determine plot size based on number of subplots, or just use use_size * rel_size
function get_plot_size(; use_size=(900,700),num_plots::Int=0, rel_size=(1,0.9))
     if num_plots > 6 # automatically determine size, this is the number of subplots
     	use_size=num_plots.*(100,100)
     elseif num_plots==2
     	use_size=(700,380)
     elseif num_plots==3
     	use_size=(700,700)
     end
     return use_size .* rel_size
end

"""
Fit f[2](ydata) as a function of f[1](xdata).

Parameters
  - `p1`, `p2`: plots to plot raw data and residuals on, respectively. Modified in place
  - `xdata`, `ydata`: x and y data points, respectively
  - `yerror`: error bars on the y data
  - `f` (kwarg): functions to apply to xdata and ydata before performing fit. Ex: To do a power law fit, take f=[log, log].
  - `xmin` (kwarg): lower cutoff for data points to include in fit
  - `cond` (kwarg): if not nothing, condition on the xdata to include in fit. Otherwise, cond will be set to x->x>=xmin
  - `use_sigma::Bool` (kwarg): whether to include yerror in weighting the least square regression (default true)
  - `plot_line::Bool` (kwarg): whether to plot a line between the data points
  - `scatter::Bool` (kwarg): whether to plot data points not included in fit as x's (default true)
  - `use_ribbon::Bool` (kwarg): whether to plot yerror as ribbon instead of error bars (do this when there are a lot of data points)
  - `label` (kwarg): label on raw data

Returns
  - `fits`: best fit parameters
  - `errs`: statistical uncertainties on best fit
"""
function plot_fit_linear(p1, p2, xdata, ydata, yerror; f=[identity,identity], xmin=1,cond=nothing, color=1, use_sigma::Bool=true, plot_line::Bool=true,
    scatter::Bool=true, use_ribbon::Bool=false, label="")

    if isnothing(cond)
        cond = x->x>=xmin
    end
    # split into data that will be used in fit, and will be excluded
    use_x = findall(cond, xdata)
    nouse_x = setdiff([1:length(xdata);], use_x)
    # get the error bars for use in fit
    if f[2]==log
        fit_errs = yerror[use_x]./ydata[use_x]
    else
	if f[2]!=identity
	    if use_sigma
	        throw(ArgumentError("Not implemented for $(f[2])"))
	    else
	        fit_errs = zeros(length(use_x))
	    end
	else
            fit_errs = yerror[use_x]
	end
    end
	    
    # plot line connecting data points, or just the data points
    if plot_line
        if yerror==zeros(length(xdata))
	   Plots.plot!(p1, xdata,ydata,color=color, label=label, msc=color)
	elseif use_ribbon
	    Plots.plot!(p1, xdata, ydata, ribbon=yerror, color=color, label=label, msc=color)
	else
           Plots.plot!(p1, xdata, ydata, yerror=yerror, color=color, label=label, msc=color)
	end
    else
        Plots.scatter!(p1, xdata[use_x], ydata[use_x], yerror=yerror[use_x], color=color, msc=color, label=label)
    end
    if scatter
        Plots.scatter!(p1, xdata[nouse_x], ydata[nouse_x], yerror=yerror[nouse_x],color=color, marker=:x, msc=color, label=label)
    end
    fits, errs, resids = fit_linear(f[1].(xdata[use_x]), f[2].(ydata[use_x]), fit_errs; use_sigma=use_sigma)
    if use_ribbon
        Plots.plot!(p2, xdata[use_x], -resids, ribbon=fit_errs,color=color, label=label)
    else
        Plots.plot!(p2, xdata[use_x], -resids, yerror=fit_errs,color=color, label=label, msc=color)
    end
    if f[2]==log && minimum(ydata .- yerror) <= 0
        Plots.ylims!(p1, (max(minimum(ydata)-0.001, 1e-5), ylims(p1)[2]))
    end
    fits, errs
end
end # module