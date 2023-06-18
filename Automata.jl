#!/usr/licensed/julia/1.7/bin/julia

"""
Express an STTI circuit as a Clifford quantum cellular automaton (CQCA).
Compute properties like recurrence times, traces over time, image of local operators, etc.
"""
module Automata

using Nemo
using Utilities
using FieldConversions
using Base.Threads
using QuantumClifford
using JLD2

DATA_DIR_SC = "/scratch/gpfs/gsommers/huse/mpt/data/"

#= Functions for writing and reading CQCA =#
export make_automaton, make_automata, read_automata

"""
Read in arrays of CQCA of the given dimensions and convert to elements of M_2a[F_2[u,u^-1]]
*Note: only implemented for qubits.*
Parameters
  - `dimensions`: array of NamedTuples (time = T, space = a) specifying the dimensions of the unit cell
  - `dir::String` (kwarg): directory to read in from
  - `automata::Dict` (kwarg): dictionary of CQCA indexed by dimensions (default empty)
  - `check::Bool` (kwarg): whether to check that automaton preserves symplectic form
  - `n_qub::Int` (kwarg): number of qubits involved in each gate (determines light cone velocity)
Returns
  - `automata`: (updated) dictionary of CQCA
"""
function read_automata(dimensions; dir = DATA_DIR_SC * "trans/automata/", automata::Dict=Dict(), check::Bool = true, n_qub::Int=2)
    R, x= LaurentPolynomialRing(GF(2), "x")
    for dims in dimensions
        automata[dims]=matrix_to_automaton.(load(dir * "$(dims.time)x$(dims.space).jld2", "automata"))
	if check
	    symplectic_2d = make_symplectic(dims.space, R)
	    for mat in automata[dims]
	        check_symplectic(mat, symplectic_2d, 1+(dims.time-1)÷dims.space*(n_qub-1))
	    end
	end
    end
    automata
end

"""
Express unit cells of Clifford gates as CQCA in matrix form. Save to file.
*Note: only implemented for qubits and two-qubit gates, but could be generalized beyond this.*
Parameters
  - `dims`: NamedTuple (time = T, space = a) specifying the dimensions of the unit cell
  - `cliffords`: array of T x a/n_qub matrices of Clifford gates
  - `dir::String` (kwarg): directory to save to
  - `automata_dict::Dict` (kwarg): dictionary that will be updated and saved to file (default empty)
  - `check::Bool`: whether to check that automaton converts correctly, and is a valid automaton
  - `n_qub::Int`: number of qubits involved in each gate (determines the lightcone velocity) (currently only works for n_qub=2)

Returns
  - `automata`: array of CQCA as elements of ``M_{2a}[F_2[u,u^{-1}]]``
"""
function make_automata(dims, cliffords; automata_dict::Dict=Dict("automata"=>[]), n_qub::Int=2, 
        check::Bool=true, dir::String="trans/automata/")
    n_qub==2 ? nothing : throw(ArgumentError("make_automata is only implemented for 2-qubit gates, you tried $(n_qub)-site gates"))
    automata = Array{Matrix}(undef, length(cliffords))
    for idx_i=1:length(cliffords)
        automata[idx_i] = make_automaton(cliffords[idx_i], check = check)
    end

    # now put into form for saving
    automata_dict["automata"] = automaton_to_matrix.(automata, 1+(dims.time-1)÷dims.space*(n_qub-1); check=check)
    save_data(automata_dict, dir * "$(dims.time)x$(dims.space)"; print=true)
    automata
end

"""
Express unit cell of Clifford gates as a CQCA in matrix form.
*Note: only implemented for qubits and two-qubit gates*

Parameters
  - `cliffords`: T x a/n_qub matrix of Clifford gates
  - `convert::Bool` (kwarg): whether to convert automaton to MatSpaceElem (rather than Matrix)
  - `check::Bool` (kwarg): whether to check that automaton preserves symplectic form
Returns
  - `automaton`: matrix with entries in ``F_2[u,u^{-1}]``, or MatSpaceElem, which encodes action of Clifford gates within one unit cell. *Columns are ordered as images of* ``X_1,X_2,...,X_a, Z_1, Z_2,...,Z_a``, *contrary to convention in paper.*
"""
function make_automaton(cliffords; convert::Bool=false, check::Bool = true)

    # how many qubits are involved in each gate? currently only implemented for 2-qubit gates
    n_qub = size(cliffords[1,1],2)
    n_qub==2 ? nothing : throw(ArgumentError("make_automata is only implemented for 2-qubit gates, you tried $(n_qub)-site gates"))

    # dimensions of unit cell: time=T, space=a
    dims = (time = size(cliffords,1), space = size(cliffords,2) * n_qub)
    R,x = LaurentPolynomialRing(GF(2),"x") 
    automaton = zeros(R, 2*dims.space,2*dims.space)
    for i=1:size(automaton,2)
        automaton[i,i]=R(1)
    end
    for t=1:dims.time
        cliff = stab_to_gf2(tensor(cliffords[t,:]...).tab)
        output = apply_gate(cliff,t)
        for x=1:2*dims.space
            old_auto = copy(automaton[x,:])
            automaton[x,:] .= R(0)
            for i=1:2*dims.space
                automaton[x,:] .+= (old_auto[i] .* output[i,:])
            end
        end
    end

    if check
        symplectic_2d = make_symplectic(dims.space, R)
	check_symplectic(automaton, symplectic_2d, 1+(dims.time-1)÷dims.space*(n_qub-1))
    end
    if convert
        S = MatrixSpace(R, size(automaton)...)
	return S(transpose(automaton))
    else
        transpose(automaton)
    end
end

#= Get the image of X's and Z's under the clifford gate cliff (tensor product of two-qubit gates). Assumes a brickwork structure where the two-site gates act inside the unit cell for odd t, and connect adjacent cells for even t =#
function apply_gate(cliff,t)
    R,x = LaurentPolynomialRing(GF(2),"x");
    output = R.(Int.(cliff))
    a = size(cliff,1)÷2
    if t%2==0 # this is where you get u^-1 and u
        for in_p in [a-1, 2*a-1]
            output[in_p, [a,2*a]] .*= x
        end
        for in_p in [a, 2*a]
            output[in_p, [a-1,2*a-1]] .*= x^(-1)
        end
        output = output[[a;1:a-1;2*a;a+1:2*a-1],[a;1:a-1;2*a;a+1:2*a-1]]
    end
    output
end

#= Functions for manipulating CQCA (reflect, shift, check valid) =#
export canonical_automaton, shift_automaton!, reflect_automaton, get_shift_mat, check_symplectic, make_symplectic

"""
Returns the symplectic 2n x 2n matrix ``\\Lambda`` over the field `R`
"""
function make_symplectic(n::Int, R)
    S = MatrixSpace(R, 2*n, 2*n)
    mat = zero(S)
    for i=1:n
        mat[i, n+i] = R(1)
	mat[n+i, i] = R(1)
    end
    mat
end

"""
Check that automaton preserves symplectic inner product: ``M' \\Lambda M = \\Lambda`` where ``\\Lambda`` is the ``2a \times 2a`` symplectic matrix and ``M' = \\overline{M^T}``
*Only works for qubits, because it calls automaton_to_matrix.*
"""
function check_symplectic(automaton, symplectic_2d, velocity::Int)
    # convert everything to matrix in case automaton or symplectic_2d are not both MatSpaceElems
    @assert Matrix(matrix_to_automaton(-automaton_to_matrix(Matrix(transpose(automaton)),velocity))) * Matrix(symplectic_2d) * Matrix(automaton) == Matrix(symplectic_2d)
end

"""
Permute columns of CQCA to bring into "canonical form" used in paper:
image of ``X_1, Z_1, X_2, Z_2,...,X_a, Z_a.``
"""
function canonical_automaton(automaton)
    space_dim = size(automaton,2)÷2
    idxs = vcat([i:space_dim:size(automaton,2) for i=1:space_dim]...)
    automaton[idxs,idxs]
end

"""
Transform the automaton so that the unit cell which used to be from qudits 1 to a is now from qudits d+1 to a+d+1. 

Parameters
  - `automaton`: element of ``M_2a[F_q[u,u^{-1}]]``, representing the action of a CQCA on unit cell a, with columns ordered ``X_1,...,X_a, Z_1,...,Z_a``
  - `d::Int`: number of sites to shift unit cell by

Returns
  - `automaton`: modified in place
"""
function shift_automaton!(automaton, d::Int)
    R,x = LaurentPolynomialRing(base_ring(automaton[1,1]),"x");
    a = size(automaton,1)÷2
    for col=[1:a-d;a+1:2*a-d]
    	for row=[a-d+1:a; 2*a-d+1:2*a]
	    automaton[row,col] *= x
	    automaton[col,row] *= x^-1
	end
    end
    automaton[[a-d+1:a;1:a-d;2*a-d+1:2*a;a+1:2*a-d],[a-d+1:a;1:a-d;2*a-d+1:2*a;a+1:2*a-d]]
end

"""
Returns the matrix M_shift, which implements a shift by d sites in a unit cell of a sites. 

Parameters
  - `a::Int`: unit cell dimension
  - `d::Int`: number of sites (<a) to shift by
  - `q::Int` (kwarg): qudit dimension (default 2)
Returns
  - `shift_mat`: MatSpaceElem M_shift, with columns ordered ``X_1, X_2,...,X_a, Z_1,Z_2,...,Z_a.``
"""
function get_shift_mat(a::Int, d::Int; q::Int = 2)
    R,x = LaurentPolynomialRing(GF(q),"x");
    S = MatrixSpace(R, 2*a, 2*a)
    shift_mat = zero(S)
    for col=1:a
        row = (col + d - 1)%a + 1
	if row < col # shifted to "next" unit cell
	    shift_mat[row,col] = x
	    shift_mat[row + a, col + a] = x
	else
	    shift_mat[row,col]=1
	    shift_mat[row+a, col+a]=1
	end
    end
    shift_mat
end

"""
Reflect the automaton about the center of the unit cell. 

Parameters
  - `automaton`: matrix representation of the CQCA, with columns ordered ``X_1,...,X_a, Z_1,...,Z_a``
  - `velocity::Int`: lightcone velocity, in units of a

Returns reflected automaton.
"""
function reflect_automaton(automaton, velocity::Int)
    space_dim = size(automaton,1)÷2
    flipped_i = [space_dim:-1:1;2*space_dim:-1:space_dim+1]
    shuffled_auto = automaton[flipped_i, flipped_i]
    return matrix_to_automaton(-automaton_to_matrix(Matrix(shuffled_auto), velocity))
end

#= Functions for getting and checking recurrence times =#
export get_periods, get_automaton_periods, has_recurrence

"""
Check whether the matrix U is equal to identity up to a shift on a system of m unit cells with PBCs

Parameters:
  - `U`: element of ``M_{2a}[F_q[u,u^{-1}]]``
  - `m::Int`: number of unit cells 
Returns: true or false whether ``U = x^d I mod (x^m+1)``
"""
function has_recurrence(U, m::Int)
    R, x= LaurentPolynomialRing(base_ring(U[1,1]), "x")
    RR = ResidueRing(R, x^m+1)
    has_recurrence(U,RR)
end

"""
Check whether the matrix U is equal to identity up to a shift in the ring RR
Parameters:
  - `U`: element of ``M_{2a}[F_q[u,u^{-1}]]``
  - `RR::Nemo.Ring`
Returns: true or false whether``U = x^d I`` (in ring RR)
"""
function has_recurrence(U,RR::Nemo.Ring)
    inv_first = 1
    try
       inv_first = inv(RR(U[1,1]))
    catch
       # not invertible in first diagonal, so not equal to shift up to identity
       return false
    end
    all(iszero(RR.(inv_first .* U .- identity_matrix(parent(U[1,1]), size(U,1)))))
end

"""
Get the period for the unitary represented by a given automaton to multiply to identity up to shifts with PBCs on systems of given lengths. 
Parameters:
  - `automaton`: element of ``M_{2a}[F_q[u,u^{-1}]]``
  - `num_cells`: system sizes (number of unit cells)
  - `max_t::Int` (kwarg): cut off to stop trying to find recurrence (default 1000)
  - `start_t::Int` (kwarg): time to start looking for recurrences (default 1)
  - `times` (kwarg): array of periods found so far (default empty). 
Returns
  - `times`: updated array of recurrence times. ith entry is `max_t`+1 if no recurrence is found for the ith system size
"""
function get_periods(automaton, num_cells; max_t::Int=1000, start_t::Int=1, times = [])
    R, x= LaurentPolynomialRing(base_ring(automaton[1,1]), "x")
    Ry, y = PolynomialRing(R, "y")
    
    res_rings = [ResidueRing(R, x^m+1) for m in num_cells]
    if isempty(times)
        times = fill(max_t + 1, length(num_cells))
    end

    remaining = Dict()
    for (length_i, m) in enumerate(num_cells)
    	if times[length_i] >= start_t # haven't found recurrence for this m yet
            remaining[m]=length_i
	    # now reset to be max_t + 1, in case I don't find recurrence
	    times[length_i] = max_t + 1
	end
    end

    # initialize to automaton U after start_t steps
    t = start_t
    U = automaton^start_t
    # store up to degree char_poly at a time
    char_poly = charpoly(Ry, automaton)
    powers = [U*automaton^(i-1) for i=1:degree(char_poly)]
    idx = 1

    # get successive powers of the automaton using recursion relation from char_poly (slightly faster than just successively multiplying?) 
    while !isempty(remaining)
        for (m, length_i) in collect(remaining)
	    if has_recurrence(powers[idx], res_rings[length_i])
                times[length_i] = t
                delete!(remaining, m)
            end
        end
	t >= max_t && break
	powers[idx] = -sum([coeff(char_poly, i) * powers[(idx+i-1)%length(powers)+1] for i=0:length(powers)-1])
        t += 1
	idx = idx%length(powers)+1
    end
    return times
end

"""
Get the recurrence times for automata on given length systems with PBCs, with unit cell a=dims.space. Return as dictionary. Save to file.
Parameters:
  - `automata`: array of CQCA all with the same unit cell dimension a, i.e. elements of ``M_{2a}[F_q[u,u^{-1}]]``
  - `lengths`: system sizes (number of qudits)
  - `dims::NamedTuple`: dimensions T, a of unit cell
  - `idxs::Array`: which of the automata to find recurrences for
  - `max_t::Int` (kwarg): cut off to stop trying to find recurrence (default 1000)
  - `start_t::Int` (kwarg): time to start looking for recurrences (default 1)
  - `periods::Dict` (kwarg): dictionary of recurrence times, indexed by automaton index 
  - `dir::String` (kwarg): directory to save data to
Returns
  - `periods`: updated dictionary of recurrence times
"""
function get_automaton_periods(automata, lengths, dims::NamedTuple, idxs; periods::Dict=Dict(),
        dir::String = "periods/", max_t::Int = 1000, start_t::Int = 1)
    period_arr = fill(max_t + 1, length(lengths), length(idxs))
    @threads for idx_i=1:length(idxs)
    	idx = idxs[idx_i]
    	# fill old data, if I have it
    	if haskey(periods, idx) && all(L->haskey(periods[idx], L), lengths)
	    period_arr[:, idx_i] = [periods[idx][L] for L in lengths]
	end
        period_arr[:,idx_i] = get_periods(automata[idx], lengths .÷dims.space,
	    max_t = max_t, start_t = start_t, times = period_arr[:,idx_i])
    end
    # now turn into dictionary (must do this outside @threads)
    for (idx_i, idx) in enumerate(idxs)
        if !haskey(periods, idx)
            periods[idx]=Dict()
        end
        for (length_i, L) in enumerate(lengths)
            periods[idx][L]=period_arr[length_i,idx_i]
        end
    end
    save_data(Dict("periods"=>periods, "tmax"=>max_t), dir * "$(dims.time)x$(dims.space)", print=true)
    periods
end

#= Functions for tracking powers of automaton or trace =#
export get_powers, get_traces, get_automaton_traces, get_pauli_image

"""
Multiply out the automaton up to time tf.
"""
function get_powers(automaton; tf=1000)
    return [automaton^t for t=0:tf]
end

"""
Get the trace of M^t, where M has characteristic polynomial char_poly, up to t=tf, using recursion relation for trace inferred from characteristic polynomial.

Parameters:
  - `automaton`: 2a x 2a matrix over ``F_q[u,u^{-1}]``
  - `char_poly`: characteristic polynomial of the automaton
  - `tf::Int` (kwarg): time up to which to compute trace (default 1000)    
Returns:
  - `traces`: array of tr(M^t) for t=0:tf
"""
function get_traces(automaton, char_poly; tf::Int=1000)
    traces = Array{typeof(tr(automaton))}(undef, tf)
    traces[1:degree(char_poly)]=[tr(automaton^i) for i=1:degree(char_poly)]
    for t=degree(char_poly)+1:tf
        traces[t] = -sum([coeff(char_poly, i)*traces[t-degree(char_poly)+i] for i=0:degree(char_poly)-1])
    end
    traces
end

"""
Get the trace of M^t, for each M in the array automata, up to t=tf, using recursion relation for trace inferred from characteristic polynomial.

Parameters:
  - `automata`: array of 2a x 2a matrix over ``F_q[u,u^{-1}]``
  - `tf::Int` (kwarg): time up to which to compute trace (default 1000)    
Returns:
  - `q::Int` (kwarg): dimension of qudit
  - `traces`: array of tr(M^t) for t=0:tf, for each automaton in the array
"""
function get_automaton_traces(automata; tf=1000, q::Int = 2)
    R, x= LaurentPolynomialRing(GF(2), "x")
    Ry, y = PolynomialRing(R, "y")
    char_polys = [charpoly(Ry, automata[idx]) for idx=1:length(automata)]
    traces = [[] for i=1:length(automata)]
    @threads for idx=1:length(automata)
        traces[idx] = get_traces(automata[idx],char_polys[idx]; tf=tf)
    end
    Array{typeof(tr(automata[1]))}.(traces)
end

"""
Get the image of the Pauli operator corresponding to the ith column of the matrix M^t, up to time tf. Express as Stabilizer tableau where the last row corresponds to time 0, first row to time tf

Parameters
  - `arrs`: powers of the automaton M, encoding time evolution up to time length(arrs)-1. Each entry of the matrix (with columns ordered in "canonical form"), originally a Laurent polynomial f, has been converted into an integer array specifying the powers of f with nonzero coefficients.
  - `i::Int`: column of matrix M^t to track the evolution of
  - `velocity::Int`: lightcone velocity in units of a (spatial dimension of unit cell)
  - `tf` (kwarg): if not nothing, max time to get the image up to. If tf is nothing, then default to tracking image up to length(arrs)-1.
  - `stabilize::Bool` (kwarg): whether to convert the final tableau into a QuantumClifford Stabilizer (default false)
  - `padding::Int` (kwarg): how much to pad the tableau by in the spatial direction (default = 0). Use padding > 0 if you want to plot the operator spreading of an operator initially supported on >1 unit cell.

Returns either a boolean matrix, or Stabilizer, giving the image of the Pauli in spacetime
"""
function get_pauli_image(arrs, i::Int, velocity; tf=nothing, stabilize::Bool=true, padding::Int=0)
    if isnothing(tf)
        tf = length(arrs)-1
    end
    space_dim = size(arrs[1,1],2)÷2
    tab = zeros(Bool, tf+1,Int(2*space_dim*(2*velocity*tf+1+2*padding)));
    for t=0:tf
        for j=1:space_dim
            # X stabilizers
            tab[end-t,(arrs[t+1][j,i] .+ Int(tf*velocity + padding)) .* space_dim .+ j] .= true
            # Z stabilizers
            tab[end-t,(arrs[t+1][j+space_dim,i] .+ Int(tf*velocity + padding)) .* space_dim .+ j .+ space_dim*Int(2*velocity*tf+1+2*padding)] .= true
        end
    end
    if stabilize
        return Stabilizer(tab)
    else
        return tab
    end
end

end # module