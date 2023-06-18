#!/usr/licensed/julia/1.7/bin/julia
"""
Conversions between finite fields, Laurent polynomials,etc.
"""
module FieldConversions

using Nemo

export matrix_to_automaton, automaton_to_matrix, poly_to_bool, from_poly, bool_to_poly, int_to_poly, to_trivariate_poly, from_trivariate_poly, ff_to_bool, bool_to_ff, nemo_to_bool, bool_to_nemo, laurent_to_array, array_to_laurent

#= Conversions between finite field (prime power) and array/matrix of bools =#
function ff_to_bool(el::fq_nmod)
    return [Bool(coeff(el, i)) for i=0:degree(parent(el))-1]
end

function ff_to_bool(V::fq_nmod_mat)
    bool_mat = [zeros(Bool, size(V,1), degree(parent(V[1,1]))) for i=1:size(V,2)]
    for i=1:size(V,1), j=1:size(V,2)
        bool_mat[j][i,:]=ff_to_bool(V[i,j])
    end
    bool_mat
end
	
function bool_to_ff(bools)
    R, z = FiniteField(2, size(bools[1], 2), "z")
    S = MatrixSpace(R, size(bools[1],1), length(bools))
    A = zero(S)
    for k=1:length(bools)
    	A[:, k] = [sum(bools[k][i,:].*[z^j for j=0:degree(R)-1]) for i=1:size(bools[k], 1)]
    end
    A
end

#= Conversion between nemo GF(2) matrix and Matrix{Bool} =#
function bool_to_nemo(V)
    @assert typeof(V[1,1])==Bool
    matrix(Nemo.GF(2), V)
end

function nemo_to_bool(V::gfp_mat)
    A = zeros(Bool, size(V)...)
    for i=1:size(V,2), j=1:size(V,1)
        A[j,i] = (V[j,i]==1)
    end
    A
end

#= Conversion between array of bools and polynomial =#

# conversion to polynomials
function bool_to_poly(bits)
    R, x = PolynomialRing(Nemo.GF(2), "x")
    R(to_poly(bits, x))
end

function int_to_poly(coefs)
    R, x = PolynomialRing(Nemo.ZZ, "x")
    R(to_poly(coefs, x))
end

# x is a variable, coefs could be integers, bools, etc
function to_poly(coefs, x)
    if isempty(coefs)
        return 0
    else
        return sum([x^(j-1)*coefs[j] for j=1:length(coefs)])
    end
end

# conversion from polynomials
function poly_to_bool(poly; convert::Bool = true, max_deg = nothing)
    coeffs = collect(coefficients(poly))
    if !isnothing(max_deg) && max_deg < length(coeffs)
        coeffs = coeffs[1:max_deg + 1]
    end
    if parent(poly).base_ring==QQ && all(isone, denominator.(coeffs))
        coeffs = numerator.(coeffs)
    end
    if !convert
        return coeffs
    elseif parent(poly).base_ring==GF(2)
        return Bool.(coeffs)
    else
        return Int.(coeffs)
    end
end

from_poly(poly; max_deg = nothing, convert::Bool = true) = poly_to_bool(poly; convert=convert,max_deg = max_deg)

"""
Conversion from integer array matrix to matrix of Laurent polynomials.
Input: `arrays`: matrix of arrays where each array corresponds to a Laurent polynomial
Output: MatSpaceElem where each entry is a Laurent polynomial
"""
function matrix_to_automaton(arrays::Matrix)
    R, x = LaurentPolynomialRing(GF(2), "x")
    S = MatrixSpace(R, size(arrays)...)
    automaton = zero(S)
    for i=1:size(arrays,2),j=1:size(arrays,1)
        automaton[j,i]=array_to_laurent(arrays[j,i])
    end
    automaton
end

"""
Conversion to integer array matrix from a matrix of Laurent polynomials.
Parameters
  - `automaton`: matrix where each entry is a Laurent polynomial
  - `max_deg`: maximum degree (negative or positive) of Laurent polynomials in the matrix
  - `check::Bool` (kwarg): if `true`, check that conversion was successful by converting back to Laurent polynomial and checking that matrix matches `automaton`. If it doesn't, you chose too small of a `max_deg.`

Returns: matrix of arrays where each array corresponds to a Laurent polynomial
"""
function automaton_to_matrix(automaton, max_deg; check::Bool=true)
    arrays = Matrix{Array{Int}}(undef, size(automaton)...)
    for i=1:size(arrays,2),j=1:size(arrays,1)
        arrays[j,i]=laurent_to_array(automaton[j,i], max_deg; check=check)
    end
    arrays
end

"""
Convert Laurent polynomial to an array of integers corresponding to the powers with nonzero coefficients.
*Only works for Laurent polynomial over ``GF(2)``.
Parameters
  - `laurent`: Laurent polynomial over ``GF(2)``
  - `max_deg`: maximum degree (negative or positive) of the Laurent polynomial
  - `check::Bool` (kwarg): if `true`, check that conversion was successful by converting back to Laurent polynomial and checking that matrix matches `automaton`. If it doesn't, you chose too small of a `max_deg.`

Returns array of powers of ``x`` that appear in the polynomial.
"""
function laurent_to_array(laurent, max_deg; check::Bool=true)
    @assert base_ring(laurent)==GF(2)
    coeffs=(-max_deg:max_deg)[findall(deg->coeff(laurent,deg)==1, -max_deg:max_deg)]
    if check
        laurent==array_to_laurent(coeffs) ? nothing : throw(
            AssertionError("max degree was chosen as $(max_deg), but $(laurent) " *
                "does not match $(array_to_laurent(coeffs))"))
    end
    coeffs
end

"""
Convert integer array `arr` into a Laurent polynomial ``\\in F_2[u,u^{-1}]``, where the term ``u^j`` appears if `j` appears in `arr`
"""
function array_to_laurent(arr)
    R, x = LaurentPolynomialRing(GF(2), "x")
    if isempty(arr) # just 0
        return R(0)
    else
        return sum([x^deg for deg in arr])
    end
end

end # module