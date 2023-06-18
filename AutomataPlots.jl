#!/usr/licensed/julia/1.7/bin/julia

"""
Plotting of operator spreading in Clifford quantum cellular automata (CQCA)
"""
module AutomataPlots

using Plots
using LaTeXStrings
using Formatting
using QuantumClifford

# my code
using PlotUtils

COLORS = palette([:white, :blue, :orange, :green])
COLOR_ARR = [:black, :blue, :orange, :green]

export plot_together, plot_separate, plot_paulis, fit_paulis, plot_trace, prep_tab, plot_multicell

"""
Get the spacetime image of the Pauli corresponding to bins, given the image of X_1,X_2,...,X_a,Z_1,...,Z_a.

Parameters
  - `bins`: length 2a boolean array, to be interpreted as a Pauli operator P on a qubits. Ex: [0, 0, 1, 0] is Z_1 I_2
  - `stabs`: length 2a array of tableaux. The first a tableaux are the spacetime images of X_1,...,X_a. The last a are images of Z_1,...,Z_a.
  - `print::Bool` (kwarg): whether to print out the string label of the Pauli operator (default true)

Returns
  - `use_tab`: spacetime image of P at time t=0,...,tf expressed as an integer matrix, where the (i,j) entry corresponds to spacetime location t=i-1, x=xdata[j], and 0=I, 1=X, 2=Z, 3=Y
  - `tag`: string label for the initial Pauli operator
"""
function prep_tab(bins, stabs; print::Bool=true)
    tag = replace("$(Stabilizer(reshape(bins, (1,length(bins)))))"[3:end], "_"=>"I")
    print && println(tag)
    s = Bool.(mod.(sum([bins[j] * stabs[j] for j=1:length(bins)]),2))
    use_tab = s[end:-1:1,1:end÷2] + s[end:-1:1,end÷2+1:end]*2
    return use_tab, tag
end

"""
Plot the operator spreading for the operator given by bins as a heatmap.
Parameters
  - `bins`: length 2a boolean array, to be interpreted as a Pauli operator P on a qubits. Ex: [0, 0, 1, 0] is ``Z_1 I_2``
  - `stabs`: length 2a array of tableaux. The first a tableaux are the spacetime images of ``X_1,...,X_a``. The last a are images of ``Z_1,...,Z_a``.
  - `space_dim::Int`: spatial dimension of unit cell.
  - `aspect` (kwarg): aspect ratio of the plot (default 1)
  - `print::Bool` (kwarg): whether to print out the string label of the Pauli operator (default `false`)
  - `colors` (kwarg): color palette to use for I, X, Z, Y (default `COLORS`)

Returns: `NamedTuple` (`tag, plt`) where
  - `tag`: string label for the initial Pauli operator
  - `plt`: heatmap of spacetime image of the given Pauli operator
"""
function plot_together(bins, stabs, space_dim; print::Bool=false, colors = COLORS, aspect = 1)
    use_tab, tag = prep_tab(bins, stabs; print=print)
    xdata = (1:size(use_tab,2)) .- (size(use_tab,2)+1-space_dim)÷2
    first_nonI = xdata[findfirst(!iszero, use_tab[1,:])]
    first_nonI == findfirst(!isequal("I"), split(tag,"")) ? nothing : throw(AssertionError(
            "$(tag), $(findfirst(!isequal("I"), split(tag, ""))), $(first_nonI)"))
    p = Plots.heatmap(xdata, 0:size(use_tab,1)-1, use_tab, aspect_ratio = aspect,
                    seriescolor=colors, clims = (0,3),colorbar=false)
    return (tag = tag, plt = p)
end

"""
Plot the operator spreading for the operator given by bins as a heatmap, with a separate panel for each site 1,...,a of the unit cell.
Parameters
  - `bins`: length 2a boolean array, to be interpreted as a Pauli operator P on a qubits. Ex: [0, 0, 1, 0] is Z_1 I_2
  - `stabs`: length 2a array of tableaux. The first a tableaux are the spacetime images of X_1,...,X_a. The last a are images of Z_1,...,Z_a.
  - `space_dim::Int`: spatial dimension of unit cell.
  - `print::Bool` (kwarg): whether to print out the string label of the Pauli operator (default false)
  - `colors` (kwarg): color palette to use for I, X, Z, Y (default COLORS)
  - `aspect` (kwarg): aspect ratio of the plot (default 1)

Returns: NamedTuple (tag, all_plots) where
  - `tag`: string label for the initial Pauli operator
  - `all_plots`: heatmaps of spacetime image of the given Pauli operator
"""
function plot_separate(bins, stabs, space_dim; print::Bool=true, colors = COLORS, aspect = 1)
    tab, tag = prep_tab(bins, stabs, print=print)
    all_plots = []
    for j=1:space_dim
        use_tab = tab[:,j:space_dim:end]
        xdata = (1:size(use_tab,2)) .- (size(use_tab,2)+1)÷2

	# print where the first identity site is at time 0
        print && println("$(j), $(use_tab[1,findfirst(isequal(0),xdata)])"); flush(stdout)
        p = Plots.heatmap(xdata, 0:size(use_tab,1)-1, use_tab, aspect_ratio = 1,
                    seriescolor=colors, clims = (0,3),colorbar=false)
        push!(all_plots, p)
    end
    return (tag = tag, plots = all_plots)
end

"""
Plot the operator spreading for the operator given by bins, with a different panel for the spacetime locations that are I, X, Y, Z
Parameters
  - `bins`: length 2a boolean array, to be interpreted as a Pauli operator P on a qubits. Ex: [0, 0, 1, 0] is Z_1 I_2
  - `stabs`: length 2a array of tableaux. The first a tableaux are the spacetime images of X_1,...,X_a. The last a are images of Z_1,...,Z_a.
  - `space_dim::Int`: spatial dimension of unit cell.
  - `velocity::Int`: light cone velocity. Also sets aspect ratio of the plot.
  - `print::Bool` (kwarg): whether to print out the string label of the Pauli operator (default false)
  - `colors` (kwarg): array of colors to use for I, X, Z, Y (default COLOR_ARR)
  - `ylabel` (kwarg): label for y axis on each panel
  - `exclude::Bool` (kwarg): whether to exclude panels where the given Pauli does not appear
  
Returns: NamedTuple (tag, all_plots) where
  - `tag`: string label for the initial Pauli operator
  - `all_plots`: heatmap showing where I, X, Y, Z are in the spacetime image of P, within lightcone
"""
function plot_paulis(bins, stabs, space_dim::Int, velocity::Int; print::Bool=false, ylabel = L"t", color_arr = COLOR_ARR,
    exclude::Bool=false)
    use_tab, tag = prep_tab(bins, stabs; print=print)
    tf = size(use_tab,1)-1
    xdata = (1:size(use_tab,2)) .- (size(use_tab,2)+1-space_dim)÷2

    # sanity checks that xdata is centered correctly
    edges = get_pauli_endpoints(use_tab[1,:], xdata, tag)

    all_plots = []
    for (pauli_i, pauli) in zip([1,2,4,3],[L"I", L"X", L"Y", L"Z"])
        tab = zeros(Bool, size(use_tab)...)
        p = Plots.plot()
        tab[findall(isequal(pauli_i-1), use_tab)] .= true
        for t=1:size(tab,1)
	    if pauli_i>1 # not identity
                @assert all(iszero,tab[t,findall(x->x<edges[1]-velocity*(t-1) || x > edges[2] + velocity*(t-1),
                            xdata)])
            else # make zero outside lightcone
                tab[t,findall(x->x<edges[1]-velocity*(t-1) || x > edges[2] + velocity*(t-1),
                        xdata)].=false
            end
        end
	if all(iszero, tab) && exclude
	    continue
	end
        Plots.annotate!((-tf*space_dim*0.6, tf/4,pauli), annotationhalign=:left,annotationfontsize=12)

        Plots.heatmap!(xdata, 0:tf, tab, aspect_ratio = velocity,
                    seriescolor=palette([:white, color_arr[pauli_i]]), 
            clims = (0,1),colorbar=false)

	# set up for a 2x2 plot
        if pauli_i in [2,3]
            Plots.plot!(p,yformatter=_->"")
        else
            Plots.plot!(p, ylabel = ylabel)
        end
        if pauli_i in [4,3]
            Plots.plot!(p, xlabel = L"x")
        else
            Plots.plot!(p, xformatter=_->"")
        end
        push!(all_plots, p)
    end
    return (tag = tag, plots = all_plots)
end

"""
Fit the fractal dimension of the I, X, Y, Z spacetime locations in the image (within the lightcone) of the spreading operator corresponding to bins.

Parameters
  - `bins`: length 2a boolean array, to be interpreted as a Pauli operator P on a qubits. Ex: [0, 0, 1, 0] is Z_1 I_2
  - `stabs`: length 2a array of tableaux. The first a tableaux are the spacetime images of X_1,...,X_a. The last a are images of Z_1,...,Z_a.
  - `space_dim::Int`: spatial dimension of unit cell.
  - `velocity::Int`: light cone velocity
  - `print::Bool` (kwarg): whether to print out the string label of the Pauli operator (default false), and some other info
  - `color_arr` (kwarg): array of colors to use for I, X, Z, Y (default COLOR_ARR)
  - `ylabel` (kwarg): label for y axis (default L"t")
  - `xmin` (kwarg): lower cutoff on time to use for fit

Returns: NamedTuple (tag, plt) where
  - `tag`: string label for the initial Pauli operator
  - `plots`: plot of the cumulative number of I,X,Y,Z inside the lightcone, with power law fits, and plot of residuals
  - `fits::Dict`: best fit parameters and errors, indexed by the label L"I", L"X", L"Y", L"Z"
  - `counts::Array`: array of counts (noncumulative) of I, X, Y, Z in each time slice
"""
function fit_paulis(bins, stabs, space_dim::Int, velocity::Int; print::Bool=false, ylabel = L"t", xmin=2^4, color_arr = COLOR_ARR)
    tabs = Dict()
    edges = []
    xdata = Dict()
    use_tab, tag = prep_tab(bins, stabs; print=print)
    xdata = (1:size(use_tab,2)) .- (size(use_tab,2)+1-space_dim)÷2
    edges = get_pauli_endpoints(use_tab[1,:], xdata, tag)
    tf = size(use_tab,1)-1    
    all_counts = []
    all_fits = Dict()
    p = Plots.plot(xscale=:log2, yscale=:log2,xlabel = L"t", ylabel = L"\sum N_{\sigma}(t')")
    p2 = Plots.plot()
    for (pauli_i, pauli) in zip([1,2,4,3],[L"I", L"X", L"Y", L"Z"])
        tab = zeros(Bool, size(use_tab)...)
        tab[findall(isequal(pauli_i-1), use_tab)] .= true
        
        for t=1:size(tab,1) # only count inside the light cone
            if pauli_i>1
                @assert all(iszero,tab[t,findall(x->x<edges[1]-velocity*(t-1) || x > edges[2] + velocity*(t-1),
                            xdata)])
            else
                tab[t,findall(x->x<edges[1]-velocity*(t-1) || x > edges[2] + velocity*(t-1),
                        xdata)].=false
            end
        end
        counts = [count(tab[t,:]) for t=1:tf+1]
	push!(all_counts, counts)
        use_t = findfirst(c->c>0, cumsum(counts)[2:end])
        if isnothing(use_t) # never appears, so don't bother plotting
            continue
        end
        xdat = (1:tf)[use_t:end]
        fits, errs = plot_fit_linear(p, p2, xdat, cumsum(counts)[use_t+1:end], zeros(tf-use_t+1); 
            f=[log,log], xmin=xmin, use_sigma = false, scatter = false, use_ribbon=true,
            color = color_arr[pauli_i], label = pauli)
        print && println(LaTeXString(format(L"d_f={:.3f}\pm{:.1e}", fits[1], errs[1])))
	all_fits[pauli] = [fits, errs]
    end
    return (tag = tag, plots = [p, p2], fits = all_fits, counts=all_counts)
end

"""
Plot the spacetime locations of nonzero coefficients in tr(M^t), up to t_f. Also fit the fractal dimension of cumulative number of nonzero coefficients.

Parameters
  - `laurent_arr`: the jth element is an array of the powers of x with nonzero coefficients in tr(M^j)
  - `tf::Int`: time up to which to plot the spacetime locations. (May be shorter than the depth to which the fit is performed.)
  - `xmin` (kwarg): lower cutoff on time to use for fit
  - `start::Int` (kwarg): time at which to start the cumulative count of nonzero coefficients (default 1)
  - `sz` (kwarg): marker size to use for scatter plot of nonzero coefficient spacetime locations. Default 0.5, but set smaller or larger depending on the depth to which you plot.
  - `ttl`: title of plot showing fit
  - `velocity::Int` (kwarg): lightcone velocity (default 1), determines the limits to use on the plot

Returns NamedTuple where
  - `plot_2d`: plot where a point at (j,t) corresponds to a nonzero coefficient in front of x^j in tr(M^t)
  - `plot_nonzero`: scatter plot of the number of nonzero coefficients as a function of time (noncumulative)
  - `fit_plots`: fit and residuals for fractal dimension
  - `fits`: best fit parameters for power law scaling of cumulative nonzero coeffs
  - `errs`: statistical errors in the best fit parameters, as returned by the fitting algorithm
"""
function plot_trace(laurent_arr, tf::Int; ttl = "", sz=0.5, start::Int=1, xmin=2^4, velocity::Int=1)

    p = Plots.plot(aspect_ratio=:equal, legend=false, xlims=(-velocity*1.02*tf,1.02*velocity*tf), ylims=(-tf/50,1.02*tf),
        title = ttl, xlabel = L"n", ylabel = L"t")
    for t=1:tf
        Plots.scatter!(laurent_arr[t], fill(t, length(laurent_arr[t])), color=:black, markersize=sz, msc=:auto)
    end
    lens = length.(laurent_arr)
    p_nonzero = Plots.scatter(lens, label="", msc=:auto, markersize=1,title = ttl,
        xlabel = L"t", ylabel = L"N(t)")
    use_start = findfirst(!iszero, lens)
    if isnothing(use_start) # trace always zero
        return (plot_2d = p, plot_nonzero = p_nonzero)
    end

    # fit to fractal dimension
    p1 = Plots.plot(xscale=:log2, yscale=:log2, xlabel = L"t", ylabel = L"\sum N(t')",
        title = ttl, legend=:topleft)
    p2 = Plots.plot(xscale=:log2, xlabel = L"t", ylabel = L"\log[\sum N(t')] -" * "fit",
        title = "Residuals",titlefontsize=12)
    start = max(start, use_start)
    weights = cumsum(lens[start:end])
    xdata = start:length(weights)+start-1
    fits, errs = plot_fit_linear(p1, p2, xdata, weights, zeros(length(weights)); 
        f=[log,log], xmin=xmin, use_sigma = false, scatter = false, use_ribbon=true)
    Plots.plot!(p1, xdata, (xdata) .^ fits[1] .* exp(fits[2]),
        label = LaTeXString(format(L"d_f={:.3f}\pm{:.1e}", fits[1], errs[1])), color=:gray)
    return (plot_2d = p, plot_nonzero = p_nonzero, fit_plots = [p1,p2], fits=fits, errs=errs)
end

"""
Plot the operator spreading for an operator that is initially supported on more than one unit cell.
Parameters
  - `binss`: array of length 2a boolean arrays, each to be interpreted as a Pauli operator P on a qubits. Length of array is number of unit cells. Ex: [[0, 0, 1, 0],[1,1,1,0]] is ``ZI(1)YX(2)``
  - `stabs`: array of length 2a arrays of tableaux. For each length 2a array, the first a tableaux are the spacetime images of ``X_1,...,X_a``. The last a are images of ``Z_1,...,Z_a``.
  - `space_dim::Int`: spatial dimension of unit cell.
  - `aspect` (kwarg): aspect ratio of the plot (default 1)
  - `print::Bool` (kwarg): whether to print out the string label of the Pauli operator (default `false`)
  - `colors` (kwarg): color palette to use for I, X, Z, Y (default `COLORS`)

Returns: `NamedTuple` (`tag, plt, tab`) where
  - `tag`: string label for the initial Pauli operator
  - `plt`: heatmap of spacetime image of the given Pauli operator
  - `tab`: 2d matrix giving the total footprint
"""
function plot_multicell(binss, stabss, space_dim; print::Bool=false, colors = COLORS, aspect = 1)
    s = Bool.(mod.(sum([sum([bins[j] * stabs[j] for j=1:length(bins)]) 
                    for (bins, stabs) in zip(binss, stabss)]),2))
    tag = prod([replace("$(Stabilizer(reshape(bins, (1,length(bins)))))"[3:end], "_"=>"I") for bins=binss])
    use_tab = s[end:-1:1,1:end÷2] + s[end:-1:1,end÷2+1:end]*2
    xdata = (1:size(use_tab,2)) .- (size(use_tab,2)+1-space_dim)÷2
    p = Plots.heatmap(xdata, 0:size(use_tab,1)-1, use_tab, aspect_ratio = aspect,
                    seriescolor=colors, clims = (0,3),colorbar=false)
    return (tag = tag, plt = p, tab = use_tab)
end

# get the left and right endpoints of the Pauli string - the first and last nonidentity sites. Check that the x coordinates are consistent.
function get_pauli_endpoints(pauli, xdata, tag)
    edges = [xdata[findfirst(!iszero, pauli)], xdata[findlast(!iszero, pauli)]]
    edges[1] == findfirst(!isequal("I"), split(tag,"")) ? nothing : throw(AssertionError(
            "Pauli $(tag) has first nonidentity at $(findfirst(!isequal("I"), split(tag, ""))), but the x coordinate is $(edges[1])"))
    edges[2] == findlast(!isequal("I"), split(tag,"")) ? nothing : throw(AssertionError(
            "Pauli $(tag) has last nonidentity at $(findlast(!isequal("I"), split(tag, ""))), but the x coordinate is $(edges[2])"))
    edges
end

end # module