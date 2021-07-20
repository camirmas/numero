### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ f5074a54-e41f-11eb-106f-4dc2c4faa110
md"""
# Chapter 9 - Gauss Elimination

## Solving for small numbers of equations

### Graphical Method

Problems when solving sets of linear equations:
- no solution (parallel)
- infinite solutions (coincident/singular)
- ill-conditioned system (lines where slopes are so close that the point of intersection is difficult to determine)

### Cramer's Rule

Uses the determinant to solve for unknowns. Singular systems have zero determinants (result is zero).

Cramer's Rule: 
> This rule states that each unknown in a system of linear algebraic equations may be expressed as a fraction of two determinants with denominator _D_ and with the numerator obtained from D by replacing the column of coefficients of the unknown in question by the constants ``b_1, b_2, ..., b_n``. For example, ``x_1`` would be computed as

```math
x_1 = \frac{\begin{equation}
   \begin{vmatrix} 
   b_1 & a_{12} & a_{13}  \\
   b_2 & a_{22} & a_{23}  \\
   b_3 & a_{32} & a_{33}  \\
   \end{vmatrix} 
\end{equation}
}{D}
```

(See Example 9.3)

### Elimination of Unknowns

To solve a pair of simultaneous equations:
1. _Forward elimination_: Manipulate the equations to eliminate one of the unknowns
2. _Back substitution_: Solve for the remaining unknown, and back-substitute the result into one of the original equations to solve for the remaining unknown

## Naive Gauss Elimination

Extends the elimination of unknowns method to large sets of equations. The naive implementation does not account for division by zero errors. 

Forward elimination:

Given the equations

```math
\begin{align}
a_{11}x_1+a_{12}x_2+a_{13}x_3+...+a_{1n}x_n = b_1 \\
a_{21}x_1+a_{22}x_2+a_{23}x_3+...+a_{2n}x_n = b_2
\end{align}
```

Multiply the first equation by ``a_{21}/a_{11}``, then subtract from the second to give

```math
(a_{22}-\frac{a_{21}}{a_{11}}a_{12})x_2+...+(a_{2n}-\frac{a_{21}}{a_{11}}a_{1n})x_n=b_2-\frac{a_{21}}{a_{11}}b_1
```

This process can be repeated for remaining equations (e.g. a third equation could be multiplied by ``a_{31}/a_{11}`` and the result subtracted from the third equation).

The first equation is called the _pivot equation_ and ``a_{11}`` is called the _pivot coefficient_.

See Example 9.5.

Observations:
1. As the system gets larger, the computation time increases greatly (roughly ``O(n^3)``)
2. Most of the effort occurs in the elimination step
"""

# ╔═╡ ff1b1f12-689c-4539-bd47-b8ce357baff1
function naive_gauss(m, b)
	rows, cols = size(m)
	
	if rows != cols
		return
	end
	
	for pivot_row=1:rows-1
		for sub_row=pivot_row+1:rows
			factor = m[sub_row,pivot_row] / m[pivot_row,pivot_row]
			
			for col=pivot_row+1:rows
				m[sub_row, col] = m[sub_row, col] - factor * m[pivot_row, col]
			end
			
			b[sub_row] = b[sub_row] - factor * b[pivot_row]
		end
	end
	
	x = zeros(rows)
	x[rows] = b[rows] / m[rows, cols]
	
	for i=rows-1:-1:1
		sum = b[i]
		
		for j=i+1:rows
			sum = sum - m[i, j] * x[j]
		end
		x[i] = sum / m[i, i]
	end
	x
end

# ╔═╡ 539b11e2-b7ff-47b2-8372-619daef7d2f3
begin
	m = [3 -.1 -.2
		 .1 7 -.3
	     .3 -.2 10]
	b = [7.85, -19.3, 71.4]
end

# ╔═╡ 15669491-8efa-4b8d-841c-e6b86e660564
naive_gauss(m, b)

# ╔═╡ 9238dd23-76e7-4d36-8383-3d2f97b5ca64
md"""
## Pitfalls

### Division By Zero

```math
\begin{align}
2x_2+3x_3=8 \\
4x_1+6x_2+7x_3=-3 \\
2x_1+x_2+6x_3=5
\end{align}
```
The technique of _pivoting_ helps with this issue (9.4.2)

### Round-Off Errors

Becomes more apparent with large numbers of equations. Rough rule of thumb, round-off is significant with 100+ equations.

### Ill-Conditioned Systems

Use the determinant to see if the solutions are ill-conditioned.

```math
\begin{align}
a_{11}x_1+a_{12}x_2=b_1 \\
a_{21}x_1+a_{22}x_2=b_2 \\ \\
x_2=-\frac{a_{11}}{a_{12}}x_1+\frac{b_1}{a_{12}} \\
x_2=-\frac{a_{21}}{a_{22}}x_1+\frac{b_2}{a_{22}} \\ \\
\frac{a_{11}}{a_{12}} \approx \frac{a_{21}}{a_{22}} \\
a_{11}a_{22} \approx a_{12}a_{21} \\ \\
a_{11}a_{22} - a_{12}a_{21} \approx 0
\end{align}
```

Example:

```math
\begin{align}
3x_1+2x_2=18 \\
-x_1+2x_2=2 \\ \\
D=3(2)-2(-1)=8
\end{align}
```
"""

# ╔═╡ 2c2182bb-5651-4104-9684-c52958eb5d5b
md"""
### Singular Systems

We can check for singularity during elimination, because the determinant of a triangular matrix is the product of its diagonal elements (see Box 9.1). If this turns out to be zero, we can terminate the program.
"""

# ╔═╡ 9513291f-7943-4b50-b74f-a5e01ebe5fad
md"""
## Improvements

### Pivoting

Switch rows so that the largest element is the pivot element.

### Scaling

If coefficients differ greatly, round-off errors can be large. Scaling mitigates this issue (See example 9.10).

**For full implementation of Gauss with pivoting, see `gauss.jl`**
"""

# ╔═╡ Cell order:
# ╟─f5074a54-e41f-11eb-106f-4dc2c4faa110
# ╠═ff1b1f12-689c-4539-bd47-b8ce357baff1
# ╠═539b11e2-b7ff-47b2-8372-619daef7d2f3
# ╠═15669491-8efa-4b8d-841c-e6b86e660564
# ╟─9238dd23-76e7-4d36-8383-3d2f97b5ca64
# ╟─2c2182bb-5651-4104-9684-c52958eb5d5b
# ╟─9513291f-7943-4b50-b74f-a5e01ebe5fad
