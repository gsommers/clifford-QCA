#!/usr/licensed/julia/1.7/bin/julia

"""
General utility functions for data analysis and manipulation
"""
module Utilities

using Statistics
using JLD2
using FileIO
using LsqFit
using StatsBase # for sample

DATA_DIR_SC = "/scratch/gpfs/gsommers/huse/mpt/data/"

export get_staggered_lengths, average_staggered, distribution_staggered, deltas_staggered, cum_hist, get_slices, save_data, get_plot_idxs, fit_linear, fit_function, find_intersection, concat, get_binaries, find_crossing, get_max_slope, fit_step, get_line, replace_zero_index, get_label, get_base, sample_sites, pad_data

# sample out of the region sample_from, either fixed number or iid
function sample_sites(num_e, sample_from::Array; fixed_n::Bool = true)
    if fixed_n
        sites = sample(sample_from, Int(num_e), replace = false)
    else
        sites = sample_from[findall(el->el<=num_e/length(sample_from), rand(length(sample_from)))]
    end
    sites
end

# sample sites on a system of length L, either a fixed number or iid
function sample_sites(num_e, L::Int; fixed_n::Bool=true)
    sample_sites(num_e, [1:L;]; fixed_n = fixed_n)
end

function get_label(labels, key; default = "b")
    if typeof(labels)<:Dict
        if haskey(labels, key)
	    label = labels[key]
	else
	    label = default
	end
    else
        label = labels
    end
    label
end

lin(x,params)=(x.*params[1]).+ params[2]

lin_origin(x, params)=lin(x,[params,0])

# returns: best fit params, stderr, residuals
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

# returns: slope, intercept, stderr, residuals
function fit_linear(xdata, ydata, sigma; p0 = [1.,0.],use_sigma::Bool=true,
	 thresh=1e-6, print::Bool=false)
    fit_function(lin, xdata, ydata, sigma, p0; use_sigma=use_sigma, thresh=thresh, print=print)
end

function get_line(xdata,ydata,sigma)
    slope = (ydata[2]-ydata[1])/(xdata[2]-xdata[1])
    err = abs(sqrt((sigma[1])^2 + (sigma[2])^2)/(xdata[2]-xdata[1]))
    intercept = ydata[1]-slope*xdata[1]
    return slope, intercept, err
end

# returns: slope, stderr, residuals
function fit_proportional(xdata, ydata, sigma; p0 = [1.],use_sigma::Bool=true)
    fit_function(lin_origin, xdata, ydata, sigma, p0; use_sigma=use_sigma)
end

# concatenate arrays along the last dimension.
function concat(arr1, arr2; force_last::Bool=false)
    sz1, sz2 = size(arr1), size(arr2)
    ndims = length(sz1)
    @assert ndims==length(sz2)
    if ndims==1 # only one dimension
        return vcat(arr1, arr2)
    elseif sz1[1:end-1]==sz2[1:end-1] # share the first n-1 dimensions
        return cat(arr1, arr2, dims=ndims)
    elseif force_last
    	# only do this if one dimension doesn't match, because you need to pad the end of the array (e.g. for entropies if they had different plateau times)
	# currently only supported for 2d else I'm fucked
    	println("Forcing concatenation along last dimension.")
	@assert sz1[2:end-1]==sz2[2:end-1]
	@assert ndims==2
	cat_arr = zeros(eltype(arr1), max(sz1[1], sz2[1]),sz1[end] + sz2[end])
	cat_arr[1:sz1[1],1:sz1[2]] = arr1
	cat_arr[1:sz2[1], sz1[2]+1:end] = arr2
	# pad to end of array
	if sz1[1] < sz2[1]
	   for i=1:sz1[2]
	       cat_arr[sz1[1]+1:end, i] = fill(cat_arr[sz1[1], i], sz2[1]-sz1[1])
	   end
        else
	   for i=sz1[2]+1:sz1[2] + sz2[2]
	       cat_arr[sz2[1]+1:end, i] = fill(cat_arr[sz2[1], i], sz1[1]-sz2[1])
	   end
	end
	return cat_arr
    else
	return vcat(arr1, arr2)
    end
end

# find closest intersection between dat1 and dat2, which are NamedTuples of arrays of x and y coords
function find_intersection(dat1, dat2;xrange=[0,1], yrange=[0,1], weight_x=10)
    x1 = dat1[:x][findall(i->dat1[:y][i]>yrange[1] && dat1[:y][i] < yrange[2] && dat1[:x][i]>xrange[1] && dat1[:x][i]<xrange[2], 1:length(dat1[:y]))]
    x2 = dat2[:x][findall(i->dat2[:y][i]>yrange[1] && dat2[:y][i] < yrange[2] && dat2[:x][i]>xrange[1] && dat2[:x][i]<xrange[2], 1:length(dat2[:y]))]

    y1 = dat1[:y][findall(el->el in x1, dat1[:x])]
    y2 = dat2[:y][findall(el->el in x2, dat2[:x])]

    min_d = Inf
    x = nothing

    for i1=1:length(y1), i2=1:length(y2)
    	d2 = (y1[i1] - y2[i2])^2 + weight_x*(x1[i1] - x2[i2])^2
	if d2 < min_d
	    min_d = d2
	    x = (x1[i1] + x2[i2])/2
	end
    end
    x
end

function find_crossing(pinned_ents, start_ent, L1, L2;best_guess=0.44, r=100, zeta=2/3)
    dat1 = dropdims(mean(pinned_ents[L1][start_ent][:,2:end], dims=1),dims=1)./L1
    dat2 = dropdims(mean(pinned_ents[L2][start_ent][:,2:end], dims=1), dims=1)./L2
    t0 = findfirst(el->el<best_guess, dat1)
    best_diff = Inf
    best_s = best_guess
    best_t = t0 / L1^zeta
    for t=max(t0-r,1):min(t0+r, length(dat1))
        s = dat1[t]
        err, t2 = findmin((L2^(zeta) ./[1:length(dat2);] .- L1^(zeta)/t).^2 .+ (dat2 .- s).^2)
        if err<best_diff
            best_diff = err
            best_s = (s + dat2[t2])/2
            best_t = (t2/L2^(zeta) + t/L1^(zeta))/2
        end
    end 
    best_s, best_t
end 

function get_plot_idxs(data, subclass, key; idxs::Dict=Dict())
    if haskey(idxs, subclass)
	return idxs[subclass]
    elseif typeof(data[subclass][key])<:Array
    	return [1:length(data[subclass][key]);]
    else
	return collect(keys(data[subclass][key]))
    end
end

function save_data(data, name;print::Bool=false)
    filename = DATA_DIR_SC * name * ".jld2"
    print && println("Saving to $(filename)")
    save(filename, data)
end

# get the length of each array, up to the first element > cutoff
function get_staggered_lengths(arrays; cutoff=0)
     if cutoff>0
	 lengths = zeros(Int64, length(arrays))
	 for (i, arr) in enumerate(arrays)
	     idx = findfirst(i->arr[i]>cutoff, [1:length(arr);])
	     if idx==nothing
		 lengths[i] = length(arr)
	     else
		 lengths[i] = idx - 1
	     end
	 end
     else
	 lengths =length.(arrays)
     end
     return lengths
 end

# get the minimum, maximum below cutoff across all arrays
function get_staggered_extrema(arrays; lengths=[])
    maxim = 1
    minim = Inf
    for (arr_i, arr) in enumerate(arrays)
    	if isempty(lengths)
	    cutoff = length(arr)
	else
	    cutoff = lengths[arr_i]
	end
	if cutoff==0
	    continue
	end
	maxm = maximum(arr[1:cutoff])
	if maxm > maxim
	    maxim = maxm
	end
	minm = minimum(arr[1:cutoff])
	if minm < minim
	    minim = minm
	end
    end
    minim, maxim
end

function pad_data(arrays, len, fill_val)
    max_len = max(maximum(length.(arrays)), len)
    data = [vcat(arr, fill(fill_val, max_len-length(arr))) for arr in arrays]
end

# what I want is a function that takes allthe samples, and then averages over them in a staggered fashion
function average_staggered(arrays; cutoff=0)
     lengths = get_staggered_lengths(arrays;cutoff=cutoff)
     max_length = maximum(lengths)
     means = zeros(max_length)
     stds = zeros(max_length)
     counts = zeros(max_length)
     idxs = [1:length(arrays);]
     for i in 1:max_length
	 filter!(idx->lengths[idx]>=i, idxs)
	 means[i] = mean([arrays[idx][i] for idx in idxs])
	 stds[i] = std([arrays[idx][i] for idx in idxs])/sqrt(length(idxs))
	 counts[i] = length(idxs)
     end
     means, stds, counts
 end

function deltas_staggered(arrays; cutoff=0)
    lengths = get_staggered_lengths(arrays; cutoff=cutoff)
    minm, maxm = get_staggered_extrema(arrays; lengths=lengths)
    deltas = [[] for i=minm:maxm]
    for (arr_i, arr) in enumerate(arrays)
        for idx in 2:lengths[arr_i]
	    delta = arr[idx] - arr[idx-1]
	    push!(deltas[arr[idx-1]-minm+1], delta)
	end
    end
    deltas, minm, maxm
end

function distribution_staggered(arrays; cutoff=0, count = 0)
    idxs = [1:length(arrays);]
    lengths = get_staggered_lengths(arrays; cutoff=cutoff)
    max_length = maximum(lengths)
    if count == 0 || count > max_length
        count = max_length
    end 
    distributions = Array{Array}(undef, count)
    for t in 1:count
        filter!(idx->lengths[idx]>=t, idxs)
	distributions[t] = [arrays[idx][t] for idx in idxs]
    end
    distributions
end

function cum_hist_minus(distribution)
    data = sort(unique(distribution))
    counts = zeros(length(data)+1)
    for (el_i, el) in enumerate(data)
        counts[el_i] = length(findall(d->d<el, distribution))
    end
    counts[end]=length(distribution)
    return push!(data, data[end]), counts/counts[end]
end

function cum_hist(distribution)
    data = sort(unique(distribution))
    counts = zeros(length(data)+1)
    for (el_i, el) in enumerate(data)
        counts[el_i+1] = length(findall(d->d<=el, distribution))
    end
    return pushfirst!(data, data[1]), counts/counts[end]
end

function get_slices(t; tail=4, numsteps=0, min_num=1)
    if numsteps>0
       step = max(min_num, Int(floor(t/numsteps)))
       return step:step:t
    end
    maxm = t-tail
    max_counts = [10, 6, 5,2,Inf];
    if t > 100
        intervals = [1,5,10,25,50]
        counts = [5,3,3,2,5]
    else
        intervals = [1,2,4,5,10,20]
    	counts = [4,3,1,4,2,4]
    end
    num = 1
    slices = []
    for (interval, count, count2) in zip(intervals, counts, max_counts)
    	rest = (maxm-num)÷interval + 1
        if rest < count2
            usecount = rest
        else
            usecount = count
        end
	push!(slices, num:interval:(num + interval*(usecount-1))...)
	num += interval * usecount
	
	if num>=maxm
	   for i=max(maxm+1,num-interval+1):(maxm+tail)
	       push!(slices, i)
	   end
	   break
	end
    end
    slices
end

function get_binaries(num, places)
    Bool.(get_base(num, places, 2))
end

# convert integer into an array of digits in a given base
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

function get_max_slope(xdata, ydata, yerror; step=1)
    @assert maximum(ydata)-minimum(ydata) >= step
    step_start=1
    step_end = length(ydata)
    fits, errs, _ = fit_linear(xdata, ydata, yerror)
    for i=1:length(ydata)-1
        end_step = findfirst(el->el>=ydata[i]+step, ydata[i+1:end])
        if isnothing(end_step)
            continue
        end
        cut_fits, cut_errs, _ = fit_linear(xdata[i:end_step+i], ydata[i:end_step+i], yerror[i:end_step+i])
        if cut_fits[1] > fits[1]
            fits = cut_fits
            errs = cut_errs
            step_start = i
            step_end = end_step+i
        end
    end
    step_start, step_end, fits, errs
end

function fit_step(xdata, ydata, yerror; lims=[0.25,0.75])
    use_me = findall(!isnan, ydata)
    _, start = findmin(y->abs(y-lims[1]),ydata[use_me])
    _, stop = findmin(y->abs(y-lims[2]),ydata[use_me])
    start = use_me[start]
    stop = use_me[stop]
    if stop-start<=1
        throw(AssertionError("Data is not finely spaced enough. start: $(start), stop: $(stop), lims: $(lims), $(ydata)"))
    end
    cut_fits, cut_errs, _ = fit_linear(xdata[start:stop], ydata[start:stop], yerror[start:stop])
    # now just get the slope, no fit
    slope, intercept, err = get_line([xdata[start], xdata[stop]], [ydata[start], ydata[stop]], [yerror[start], yerror[stop]])
    return (start=start, stop=stop, fits=cut_fits, errs=cut_errs, line=[slope, intercept], slope_err = err)
end

function replace_zero_index(idx, arr, j; default=-1)
    if idx>0
        return arr[idx, j]
    elseif idx==0 # access before the start of the array: just set to default value
        return default
    else
        throw(ArgumentError("Gave negative index to array: $(idx)"))
    end
end

end