### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 9052fffd-fb04-4aa6-a476-22c753fd0c74
include("derivative.jl")

# ╔═╡ 3d81fb46-e1ad-11eb-04ec-e50f53e65178
md"# Chapter 6 - Open Methods"

# ╔═╡ 878916b8-7e02-47db-a509-cfc20b240e19
md"""
- _Open methods_ require only a single starting value of _x_ or two starting values that do not necessarily bracket the root
- They sometimes _diverge_ or move away from the true root as the computation progresses
- When they do converge, they are usually faster than bracketing methods
"""

# ╔═╡ bf709f85-c629-4f66-8e3f-85461ca344b9
md"""
## Simple Fixed-Point Iteration

Method: put _x_ on one side of the equation, or add _x_ to both sides. For example,

```math
x^2-2x+3=0
```

becomes
```math
x=\frac{x^2+3}{2}
```

whereas
```math
sin(x)=0
```

becomes
```math
x=sin(x)+x
```

Again, the approximate error is
```math
\epsilon_a = \left\lvert\frac{x_{i+1}-x_i}{x_{i+1}}\right\rvert100\%
```

**Linear convergence** - The true percent relative error is proportional to the error from the previous iteration. Convergence occurs if `|g'(x)|<1`, when the magnitude of the slope of `g(x)` is less than the slope of the line `f(x) = x`.
"""

# ╔═╡ 46531212-8ee5-431e-9b79-9953f61de69e
function fixed_point(f, x_0, e_s, i_max)
	x_r = x_0
	e_a = 0
	
	for i=1:i_max
		x_r_old = x_r
		x_r = f(x_r_old)
		
		if x_r != 0
			e_a = abs((x_r - x_r_old) / x_r) * 100
		end
		
		if e_a < e_s
			return x_r
		end
	end
	
	return x_r
end

# ╔═╡ bd53783e-0bbf-4702-b226-a39b9827fb3c
fixed_point((x) -> exp(-x), 0, .5, 10)

# ╔═╡ 01ab626c-1e13-40ec-b4fa-b27c9a18fb69
md"""
## The Newton-Raphson Method

- Most widely used root-locating formula
- If the initial guess is ``x_i``, a tangent can be extended from the point ``[x_i, f(x_i)]``. The point where this tangent crosses the x axis gives an improved root estimate.

The first derivative at ``x`` is equivalent to the slope:

```math
f'(x_i)=\frac{f(x_i)-0}{x_i-x_{i+1}}
```

which rearranges to

```math
x_{i+1}=x_i-\frac{f(x_i)}{f'(x_i)}
```

This can also be derived using a Taylor Series (see Box 6.2).

Pitfalls:
- Slow for functions that have an inflection point (``f''(x)=0``) near a root
- Multiple roots can cause the function to jump to farther away roots
- Zero slopes cause division by zero errors
"""

# ╔═╡ 523abe51-fc23-44ac-8d51-91416a1ce8b1
newton_raphson_guess(f, x_i) = x_i - f(x_i)/derivative(f, x_i)

# ╔═╡ ccafedae-ae85-41ee-929d-3fb907f28540
function newton_raphson(f, x_0, e_s, i_max)
	x_r = x_0
	e_a = 0
	
	for i=1:i_max
		x_r_old = x_r
		x_r = newton_raphson_guess(f, x_r_old)
		
		if x_r != 0
			e_a = abs((x_r - x_r_old) / x_r) * 100
		end
		
		if e_a < e_s
			return x_r
		end
	end
	
	return x_r
end

# ╔═╡ 8e186197-e9fd-4562-bfb5-6eecc4a4be2d
newton_raphson((x) -> exp(-x)-x, 0, .5, 10)

# ╔═╡ c38234d7-cd49-4897-94c3-d4199b7a9110
md"""
## The Secant Method

Uses a backward finite divided difference to solve for difficult derivatives. Requires two initial guesses, similar to bracketing. Similar to the false-position method, but it can diverge.

```math
f'(x_i) \approx \frac{f(x_{i-1})-f(x_i)}{x_{i-1}-x_i}
```

Substituted into Newton-Raphson,

```math
x_{i+1} = x_i - \frac{f(x_i)(x_{i-1}-x_i)}{f(x_{i-1})-f(x_i)}
```
"""

# ╔═╡ c91de686-42f1-42a2-9793-b533e8fb62bc
md"""
## Brent's Method

Uses open method where possible, but reverts to bracketing if necessary. Uses bisection for bracketing. Uses secant method and inverse quadratic interpolation for open methods.

#### Inverse Quadratic Interpolation

Similar to secant, but uses a quadratic function that goes through three points. It is "inverse" because it is defined in terms of ``y`` (``x=f(y)``), avoiding issues where the parabola may have complex roots.

```math
x_{i+1} = \frac{y_{i-1}y_i}{(y_{i-2}-y_{i-1})(y_{i-2}-y_i)}x_{i-2} + \frac{y_{i-2}y_i}{(y_{i-1}-y_{i-2})(y_{i-1}-y_i)}x_{i-1} + \frac{y_{i-2}y_{i-1}}{(y_i-y_{i-2})(y_i-y_{i-1})}x_i
```

Does not work if the three ``y`` values are not distinct. In such a case, fall back to secant method. If ``y_{i-2}=y_{i-1}``, use secant method with ``x_{i-1}, x_i``. If ``y_{i-1}=y_i``, use ``x_{i-2}, x_{i-1}``.

TODO: Implementation
"""

# ╔═╡ Cell order:
# ╟─3d81fb46-e1ad-11eb-04ec-e50f53e65178
# ╟─878916b8-7e02-47db-a509-cfc20b240e19
# ╟─bf709f85-c629-4f66-8e3f-85461ca344b9
# ╠═46531212-8ee5-431e-9b79-9953f61de69e
# ╠═bd53783e-0bbf-4702-b226-a39b9827fb3c
# ╟─01ab626c-1e13-40ec-b4fa-b27c9a18fb69
# ╠═9052fffd-fb04-4aa6-a476-22c753fd0c74
# ╠═523abe51-fc23-44ac-8d51-91416a1ce8b1
# ╠═ccafedae-ae85-41ee-929d-3fb907f28540
# ╠═8e186197-e9fd-4562-bfb5-6eecc4a4be2d
# ╟─c38234d7-cd49-4897-94c3-d4199b7a9110
# ╟─c91de686-42f1-42a2-9793-b533e8fb62bc
