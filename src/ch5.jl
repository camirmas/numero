### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 9100dc36-de94-11eb-3dca-41060ad46b86
using Plots

# ╔═╡ e5e78310-2019-4070-a549-19319245d5ce
md"""
# Chapter 5 - Bracketing Methods
"""

# ╔═╡ 07a28fae-85cd-4ebd-8576-049e0191fc7f
f(x) = sin(10x) + cos(3x)

# ╔═╡ bfa8330f-a10f-4177-a8cc-14d6dfaec84e
begin
	plot(f, xlims=(0, 5))
end

# ╔═╡ 7f0c6337-5917-45b2-bdc2-12d7e2b8fc05
md"""
## The Bisection Method

1. Choose lower `x_l` and upper `x_u` guesses for the root such that the function changes sign over the interval. This can be checked by ensuring that `f(x_l)f(x_u) < 0`.
2. An estimate of the root `x_r` is determined by `x_r = (x_l+x_u)/2`
3. Make the following evaluations to determine in which subinterval the root lies:
    * If `f(x_l)f(x_r) < 0`, the root lies in the lower subinterval. Therefore, set `x_u = x_r` and return to step 2.
    * If `f(x_l)f(x_r) > 0`, the root lies in the upper subinterval. Therefore, set `x_l = x_r` and return to step 2.
    * If `f(x_l)f(x_r) = 0`, the root equals `x_r`; terminate the computation.

The computation can be ended after the approximate error falls below a particular threshold. For this method, the approximate error will always be above the actual.

```math
\epsilon_a = \left\lvert\frac{x_r^{new}-x_r^{old}}{x_r^{new}}\right\rvert100\%
```
"""

# ╔═╡ 9513dc83-fb0d-4813-a5fe-ebf465c2a587
function bisect(f, x_l, x_u, x_r, i_max, e_s)
	if f(x_l)*f(x_u) >= 0
		return
	end

	f_l = f(x_l)
		
	for i=0:i_max
		x_r_old = x_r
		x_r = (x_l + x_u)/2
		f_r = f(x_r)
		
		if x_r != 0
			e_a = abs((x_r - x_r_old) / x_r) * 100
		end
		
		test = f_l * f_r

		if test < 0
			x_u = x_r
		elseif test > 0
			x_l = x_r
			f_l = f_r
		else
			e_a = 0
		end
		
		if e_a < e_s
			return x_r
		end
	end
	
	return x_r
end

# ╔═╡ 352dcd87-b2dc-4b12-b33d-c111d9609370
fn(c) = 668.06/c*(1-exp(-.146843c))-40

# ╔═╡ 45509df2-f2ef-48ec-b9f2-10dadc12be52
fn(14.75)

# ╔═╡ 6955eb8e-6bc9-4263-93bf-1abff129594d
bisect(f, 12, 16, 14.75, 10, .5) 

# ╔═╡ cf7baab5-5692-45b1-a67d-4a39ceb95be2
md"""
## The False-Position Method

An alternative method based on graphical insight. Join f(x_l) and f(x_u) by a straight line. The intersection of this line with the x axis gives an improved estimate of the root.

This method is often more effective than bisection, but not always. In addition to using the approximate error calculations, the results should always be checked by substituting the root estimate back into the original equation and determining if the result is close to zero.


```math
x_r = x_u - \frac{f(x_u)(x_l-x_u)}{f(x_l)-f(x_u)}
```
"""

# ╔═╡ ea691cf1-d206-4d47-a5a0-43b2456f5582
false_position_guess(f, x_l, x_u) = x_r = x_u - f(x_u)*(x_l-x_u) / (f(x_l)-f(x_u))

# ╔═╡ c8d98f6b-94aa-4415-927c-1f04dd44c29d
false_position_guess(fn, 12, 16)

# ╔═╡ f9e09248-711c-46b3-a6a3-efc280b7ba35
md"""
In order to mitigate the "one-sidedness" of this method (i.e. one of the bracketing points will tend to stay fixed), the stagnant bound can be modified.
"""

# ╔═╡ 3db9e84f-2c5b-4306-9bcd-36440d2fc53a
function false_position(f, x_l, x_u, e_s, i_max, x_r)
	f_l = f(x_l)
	f_u = f(x_u)
	i_l = 0
	i_u = 0
	
	for i=1:i_max
		x_r_old = x_r
		x_r = false_position_guess(f, x_l, x_u)
		f_r = f(x_r)
		
		if x_r != 0
			e_a = abs((x_r - x_r_old) / x_r) * 100
		end
		
		test = f_l * f_r
		
		if test < 0
			x_u = x_r
			f_u = f(x_u)
			i_u = 0
			i_l += 1
			
			if i_l >= 2
				f_l /= 2
			end
		elseif test > 0
			x_l = x_r
			f_l = f(x_l)
			i_l = 0
			i_u += 1
			
			if i_u >= 2
				f_u /= 2
			end
		else
			e_a = 0
		end
		
		if e_a < e_s
			return x_r
		end
	end
	
	return x_r
end

# ╔═╡ 58faf484-9abb-4911-80c4-0ad2079c4969
res = false_position(f, 12, 16, .5, 10, 14.75)

# ╔═╡ edd4d422-33b1-4e7e-bea8-5bb21d43ea21
f(res)

# ╔═╡ 134a9e04-b998-4b8f-a853-46f8305a7513
f(res)

# ╔═╡ Cell order:
# ╟─e5e78310-2019-4070-a549-19319245d5ce
# ╠═9100dc36-de94-11eb-3dca-41060ad46b86
# ╠═07a28fae-85cd-4ebd-8576-049e0191fc7f
# ╠═bfa8330f-a10f-4177-a8cc-14d6dfaec84e
# ╟─7f0c6337-5917-45b2-bdc2-12d7e2b8fc05
# ╠═9513dc83-fb0d-4813-a5fe-ebf465c2a587
# ╠═352dcd87-b2dc-4b12-b33d-c111d9609370
# ╠═45509df2-f2ef-48ec-b9f2-10dadc12be52
# ╠═6955eb8e-6bc9-4263-93bf-1abff129594d
# ╠═edd4d422-33b1-4e7e-bea8-5bb21d43ea21
# ╟─cf7baab5-5692-45b1-a67d-4a39ceb95be2
# ╠═ea691cf1-d206-4d47-a5a0-43b2456f5582
# ╠═c8d98f6b-94aa-4415-927c-1f04dd44c29d
# ╟─f9e09248-711c-46b3-a6a3-efc280b7ba35
# ╠═3db9e84f-2c5b-4306-9bcd-36440d2fc53a
# ╠═58faf484-9abb-4911-80c4-0ad2079c4969
# ╠═134a9e04-b998-4b8f-a853-46f8305a7513
