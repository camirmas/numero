### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 9c382188-ecb1-11eb-3843-258ac1439564
md"""
# Ch 13 - One-Dimensional Unconstrained Optimization

_Optimization_ involves searching for either the minimum or maximum. The optimum is the point where the curve is flat (i.e. where ``f'(x)=0``). The second derivative, ``f''(x)``, indicates whether the optimum is a minimum (``f''(x) > 0``) or a maximum (``f''(x) < 0``).

_Multimodal_ functions have both local and global optima.

Just as with root finding, optimization in one dimension can be divided into bracketing and open methods. Open methods tend to be faster but can also diverge.

## Golden-Section Search

Similar to bisection approach to root finding. Verify that a root exists by checking that ``f(x_l)`` and ``f(x_u)`` have different signts. Use midpoint to guess, then narrow the window.

Midpoint:
```math
x_r = \frac{x_l+x_u}{2}
```

Interior points are determined using the _golden ratio_:

```math
\begin{align}
d&=\frac{\sqrt{5}-1}{2}(x_u-x_l) \\
x_1&=x_l+d \\
x_2&=x_u-d
\end{align}
```

1. If ``f(x_1) > f(x_2)``, then the domain of ``x`` to the left of ``x_2``, from ``x_l`` to ``x_2``, can be eliminated. ``x_2`` becomes the new ``x_l``.
2. If ``f(x_2) > f(x_1)``, then the domain of ``x`` to the right of ``x_1``, from ``x_1`` to ``x_u`` would be eliminated. ``x_1`` becomes the new ``x_u``.

The benefit of the golden ratio is that we do not have to recalculate all of the function values for the next iteration, as we are reusing one of them.

Estimated error:
```math
\epsilon_a=(1-R)\left\lvert\frac{x_u-x_l}{x_{opt}}\right\rvert100\%
```
"""

# ╔═╡ acd57efd-288a-4465-ab75-0c96b1fc0d4f
function gold(x_l, x_u, max_it, es, f, comp=(>))
	R=(5^(.5)-1)/2 # golden ratio
	
	d = R * (x_u-x_l)
	x_1 = x_l + d # inner interval 1
	x_2 = x_u - d # inner interval 2
	
	f_1 = f(x_1)
	f_2 = f(x_2)
	
	if comp(f_1, f_2)
		x_opt = x_1
		f_x = f_1 # value of f(x) at the current optimum
	else
		x_opt = x_2
		f_x = f_2
	end
	
	for i=1:max_it
		d = R*d
		x_int = x_u - x_l
		
		if comp(f_1, f_2)
			x_l = x_2
			x_2 = x_1
			x_1 = x_l + d
			f_2 = f_1
			f_1 = f(x_1)
		else
			x_u = x_1
			x_1 = x_2
			x_2 = x_u - d
			f_1 = f_2
			f_2 = f(x_2)
		end
		
		if comp(f_1, f_2)
			x_opt = x_1
			f_x = f_1
		end
		
		if x_opt != 0
			ea = (1-R)*abs(x_int/x_opt)*100 # estimated error
		end
		
		if ea < es
			return x_opt
		end
	end
	return x_opt
end

# ╔═╡ 37bb1c18-1b9e-4f47-b23d-8c6fd0250263
begin
	f(x) = 2sin(x)-x^2/10
	gold(0, 4, 8, .01, f)
end

# ╔═╡ 69ee8dd3-3b22-4789-b996-fc0dff0d5875
md"""
## Parabolic Interpolation

Find a new point

```math
x_3 = \frac{f(x_0)(x_1^2-x_2^2)+f(x_1)(x_2^2-x_0^2)+f(x_2)(x_0^2-x_1^2)}{2f(x_0)(x_1-x_2)+2f(x_1)(x_2-x_0)+2f(x_2)(x_0-x_1)}
```

then apply a bracketing approach similar to bisection or golden-section.
"""

# ╔═╡ 25816687-42b2-46e0-b058-bbffeb51e32f
function parabolic(x_0, x_1, x_2, f, max_it, es)
	f_0 = f(x_0)
	f_1 = f(x_1)
	f_2 = f(x_2)
	
	for _=1:max_it
		
	end
	
	x_3 = (f_0(x_1^2-x_2^2)+f_1(x_2^2-x_0^2)+f_2(x_0^2-x_1^2)) /
		(2f_0(x_1-x_2)+2f_1(x_2-x_0)+2f_2(x_0-x_1)) 
end

# ╔═╡ Cell order:
# ╠═9c382188-ecb1-11eb-3843-258ac1439564
# ╠═acd57efd-288a-4465-ab75-0c96b1fc0d4f
# ╠═37bb1c18-1b9e-4f47-b23d-8c6fd0250263
# ╠═69ee8dd3-3b22-4789-b996-fc0dff0d5875
# ╠═25816687-42b2-46e0-b058-bbffeb51e32f
