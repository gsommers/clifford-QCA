#!/usr/licensed/julia/1.7/bin/julia
"""
Label the 36 iSWAP-core circuits on the square lattice.
"""
module DualUnitary

using Formatting
using LaTeXStrings

export get_idx_1d, get_idx_2d, get_idx_4d, get_gate_label

C1_GATES = [L"I", L"R_X[\pi/2]", L"R_{-}[2\pi/3]",
        L"R_{+}[2\pi/3]", L"R_Z[\pi/2]", L"R_Y[\pi/2]"]
S1_IDXS = [1,2,6]

get_idx(i; num=6)=(i-1)%num^2+1

"""
Convert from 1d idx in a num^2 length list to a 2d idx in a num x num matrix. num defaults to 6, the number of possible one-site gates on each leg.
"""
function get_idx_2d(idx; num=6)
    i = get_idx(idx; num = num)
    return [(i - 1)÷num+1, (i - 1)%num+1]
end

"""
Convert from n-d index in a matrix where the first n-1 dimensions are given by num, to a 1d index.
""" 
function get_idx_1d(idx; num=[6])
    return sum([prod(num[i:end])*(idx[i]-1) for i=1:length(num)])+idx[end]
end

"""
Convert from a 1d index into a 4d index. Useful for labeling a gate decorated with a one-site gate on each leg.

Parameters
  - `i::Int`: 1d index
  - `nums` (kwarg): number of choices on each leg/length of each of the four dimensions. Default [6,6,3,3] for the ways to decorate an iSWAP or CNOT core.
  - `make_s1::Bool` (kwarg): if true, convert the indices on the outgoing legs (3 distinct choices) into indices into the array of all 6 gates

Returns 4d index
"""
function get_idx_4d(i::Int; nums = [6,6,3,3], make_s1::Bool=true)
    coords=[((i-1)%prod(nums[end:-1:j]))÷prod(nums[end:-1:j+1]) + 1 for j=1:length(nums)]
    if make_s1
        return [coords[1], coords[2], S1_IDXS[coords[3]], S1_IDXS[coords[4]]]
    else
        return coords
    end
end

"""
Get the label for the one site gates of the gate with index idx.
"""
function get_gate_label(idx; nums = [6,3], make_4d::Bool=false)
    if make_4d
       idx_4d = get_idx_4d(idx, (nums[1], nums[1], nums[2], nums[2]))
       return format("({}, {}, {}, {})", C1_GATES[idx_4d[1]], C1_GATES[idx_4d[2]], C1_GATES[S1_IDXS[idx_4d[3]]], C1_GATES[S1_IDXS[idx_4d[4]]])
    else
       idx_2d = get_idx_2d(idx; num=nums[1])
       return LaTeXString(format("({},{})", C1_GATES[idx_2d[1]], C1_GATES[idx_2d[2]]))
    end
end

end # module

