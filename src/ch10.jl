### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ ab0aaede-ea4e-11eb-328e-5f12ec812d85
md"""
# Chapter 10 - _LU_ Decomposition and Matrix Inversion

## _LU_ Decomposition

- Useful when many right-hand side vectors {B} must be computed for the same matrix of coefficients [A]
- Requires pivoting to avoid division by zero

Start with

``[A]\{X\}-\{B\}=0``

then use elimination to get an upper triangular system, represented as

``[U]\{X\}-\{D\}=0``

then assume there is a lower diagonal matrix with 1's on the diagonal, such that

``[L]\{[U]\{X\}-\{D\}\}=[A]\{X\}-\{B\}``

thus, by the rules of matrix multiplication,

``[L][U]=[A]``

and

``[L]\{D\}=[B]``

Strategy
1. _LU Decomposition_. [A] is factored, or "decomposed," into lower [L] and upper [U] triangular matrices
2. _Substitution_. [L] and [U] are used to determine a solution {X} for the right-hand-side {B}. First, use `[L]{D}={B}` to generate an intermediate vector {D} by forward substitution. Then, the result is substituted into `[U]{X}-{D}=0`, which can be solved by back substitution for {X}.

Gauss elimination can be used for decomposition. ``[U]`` is a direct product of forward elimination.

``[L]`` is formed using elimination factors (See. 10.1.2).
"""

# ╔═╡ edb3a016-61a2-48ab-b8bf-4312a476f129
"""
Decomposition phase (naive, does not include pivoting).
For a full pseudocode implementation, see Figure 10.2).
"""
function decompose(a)
	n, _ = size(a)
	
	for k=1:n-1
		for i=k+1:n
			factor = a[i, k]/a[k,k]
			a[i, k] = factor
			
			for j=k+1:n
				a[i,j] -= factor*a[k,j]
			end
		end
	end
	
	return a
end

# ╔═╡ 6f079425-9ede-4e49-b3b9-3fb7f5832081
a = [3 -0.1 -0.2
	0.1 7 -0.3
	0.3 -0.2 10]

# ╔═╡ 556eea55-5796-4f71-99ae-82b547696e60
decompose(a)

# ╔═╡ 7ce699d7-9f72-42db-afbd-9f565c7500fd
function substitute(a, b)
	n, _ = size(a)
	x = zeros(n)
	
	# forward
	for i=2:n
		sum = b[i]
		
		for j=1:i-1
			sum -= a[i, j]*b[j]
		end
		b[i] = sum
	end
	
	# backward
	x[n] = b[n]/a[n,n]
	
	for i=n-1:-1:1
		sum = 0
		
		for j=i+1:n
			sum += a[i,j]*x[j]
		end
		x[i] = (b[i]-sum)/a[i,i]
	end
end

# ╔═╡ 0e8fb5b9-5a94-45ce-8d46-69be3ed87dd3
md"""
## The Matrix Inverse

The inverse can be used to provide a solution:

``\{X\}=[A]^{-1}\{B\}``

It can be computed by generating solutions with unit vectors as the right-hand-side-constants. _LU_ Decomposition is the best way to implement this, as it provides an efficient means to evaluate multiple right-hand-side vectors. (See Example 10.3)

An additional property of the inverse is that it represents the response of a single part of the system to a unit stimulus of any other part of the system.

**Properties of linear systems**

_Superposition_: if a system is subject to several different stimuli (the _b_'s), the responses can be computed individually and the results summed to obtain a total response.

_Proportionality_: multiplying the stimuli by a quantity results in the response to those stimuli being multiplied by the same quantity.

## Error Analysis and System Condition

The inverse can be used to test if systems are ill-conditioned:
1. Scale the matrix of coefficients ``[A]`` so that the largest element in each row is 1. Invert the scaled matrix and see if there are elements of ``[A]^-1`` that are several orders of magnitude greater than 1. If so, this indicates ill-conditioning
2. Multiply the inverse by the original coefficient matrix and assess whether the result is close to the identity matrix. If not, this indicates ill-conditioning
3. Invert the inverted matrix and assess whether the result is sufficiently close to the original coefficient matrix. If not, this indicates ill-conditioning
"""

# ╔═╡ Cell order:
# ╟─ab0aaede-ea4e-11eb-328e-5f12ec812d85
# ╠═edb3a016-61a2-48ab-b8bf-4312a476f129
# ╠═6f079425-9ede-4e49-b3b9-3fb7f5832081
# ╠═556eea55-5796-4f71-99ae-82b547696e60
# ╠═7ce699d7-9f72-42db-afbd-9f565c7500fd
# ╠═0e8fb5b9-5a94-45ce-8d46-69be3ed87dd3
