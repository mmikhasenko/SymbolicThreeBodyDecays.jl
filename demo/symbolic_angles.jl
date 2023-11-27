### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ c6b0982a-d3ba-4736-8798-12fb497e57dc
# ╠═╡ show_logs = false
begin
    cd(joinpath(@__DIR__, ".."))
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()
    # 
    using SymbolicThreeBodyDecays
    using ThreeBodyDecay
    using SymPy
end

# ╔═╡ eec59aa0-8c9d-4ab8-b34b-e86489a6c1cc
md"""
# Angles with Symbolic expressions

The notebook calls `ThreeBodyDecay` implementation passing `SymPy` object.
See Appendix of [Dalitz-Plot Decomposition paper](https://inspirehep.net/literature/1758460) for details.
"""

# ╔═╡ abc1ba70-f1a9-4698-8a96-774d24f063ce
begin
    @syms m1 m2 m3 m0
    @syms σ1 σ2 σ3
end;

# ╔═╡ da366cd6-6204-4ed3-abb5-3db04954449a
ms = ThreeBodyMasses(m1, m2, m3; m0)

# ╔═╡ 6931fe28-767a-4bbf-99bc-2a70b83d409c
σs = Invariants(; σ1, σ2, σ3)

# ╔═╡ a511b9ff-c04b-453f-b261-24846b711679
## Formatting
begin
	latex(x::Pair{String,Sym}) = "$(x[1]) &= $(sympy.latex(x[2]))"
	latex(x::Pair{String,String}) = "$(x[1]) &= $(x[2])"
	# 
	format_equation(x::Pair) = Markdown.parse(
	"""
	```math
	\\begin{align}
	$(latex(x))
	\\end{align}
	```
	""")
	# 
	format_equation(x::Vector{T} where {T<:Pair}) = Markdown.parse(
	"""
	```math
	\\small
	\\begin{align}
	$(join(latex.(x), "\\,,\\\\"))
	\\end{align}
	```
	""")
end

# ╔═╡ a7f04b3c-7711-4111-834d-e6ac50b7a7d0
md"""
## Scattering angles
"""

# ╔═╡ 7f0ea9e8-5b8c-4c96-b958-668671aa92bf
map(1:3) do k
	i,j = ijk(k)
	"\\cos\\color{green}\\theta_{$i$j}" => cosθij(k,σs, ms^2)
end |> format_equation

# ╔═╡ 1c867f6f-842a-4ebc-8ea3-892033ffdfe5
md"""
## Wigner rotations
Green color stands of positive angles, while the negative angles are colored red.
1) All indices are different
"""

# ╔═╡ 173e97dd-f593-4db1-a927-68401a34351e
map(1:3) do k
	i,j = ijk(k)
	_wr = wr(i,j,k)
	c = ispositive(_wr) ? "green" : "red"
	"\\cos\\color{$c}\\zeta_{$i$j}^$k" => cosζ(_wr,σs, ms^2)
end |> format_equation

# ╔═╡ 9465baa1-58e0-4d0c-a4c2-3a71e6eea060
md"""
2) Two indices are the same
"""

# ╔═╡ 9c3eb17f-ee50-4daa-b6da-5676ace61cce
map(1:3) do k
	i,j = ijk(k)
	_wr = wr(i,k,k)
	c = ispositive(_wr) ? "green" : "red"
	"\\cos\\color{$c}\\zeta_{$i$k}^$k" => cosζ(_wr,σs, ms^2)
end |> format_equation

# ╔═╡ a858f640-9342-4337-a4b7-3bab6484fee4
map(1:3) do k
	i,j = ijk(k)
	_wr1 = wr(i,k,k); c1 = ispositive(_wr1) ? "green" : "red"
	_wr2 = wr(k,i,k); c2 = ispositive(_wr2) ? "green" : "red"
	@assert cosζ(_wr1,σs, ms^2) == cosζ(_wr2,σs, ms^2)
	"\\color{$c1}\\zeta_{$i$k}^$k" => "-\\color{$c2}\\zeta_{$k$i}^$k"
end |> format_equation

# ╔═╡ 33c922d5-3f11-49f7-83a7-c58800fe130f
md"""
3) Mother particle angles
"""

# ╔═╡ de63b6f5-b534-4f3e-8bd6-149bf8738ef3
map(1:3) do k
	i,j = ijk(k)
	_wr = wr(i,j,0)
	c = ispositive(_wr) ? "green" : "red"
	"\\cos\\color{$c}\\zeta_{$i$j}^0" => cosζ(_wr,σs, ms^2)
end |> format_equation

# ╔═╡ 216b21f9-38f3-4ce0-a648-8d9a1171b4b6
md"""
4) Trivial angles
"""

# ╔═╡ 5aad199f-be14-4b9d-b1f0-1e13f553f6fc
vcat(map(1:3) do k
	i,j = ijk(k)
	[
	"\\cos\\zeta_{$k$k}^$k" => cosζ(wr(k,k,k),σs, ms^2),
	"\\cos\\zeta_{$i$i}^$k" => cosζ(wr(i,i,k),σs, ms^2),
	"\\cos\\zeta_{$j$j}^$k" => cosζ(wr(j,j,k),σs, ms^2),
	"\\cos\\zeta_{$k$k}^0" => cosζ(wr(k,k,0),σs, ms^2)]
end...) |> format_equation

# ╔═╡ Cell order:
# ╟─eec59aa0-8c9d-4ab8-b34b-e86489a6c1cc
# ╠═c6b0982a-d3ba-4736-8798-12fb497e57dc
# ╠═abc1ba70-f1a9-4698-8a96-774d24f063ce
# ╠═da366cd6-6204-4ed3-abb5-3db04954449a
# ╠═6931fe28-767a-4bbf-99bc-2a70b83d409c
# ╟─a511b9ff-c04b-453f-b261-24846b711679
# ╟─a7f04b3c-7711-4111-834d-e6ac50b7a7d0
# ╟─7f0ea9e8-5b8c-4c96-b958-668671aa92bf
# ╟─1c867f6f-842a-4ebc-8ea3-892033ffdfe5
# ╟─173e97dd-f593-4db1-a927-68401a34351e
# ╟─9465baa1-58e0-4d0c-a4c2-3a71e6eea060
# ╟─9c3eb17f-ee50-4daa-b6da-5676ace61cce
# ╟─a858f640-9342-4337-a4b7-3bab6484fee4
# ╟─33c922d5-3f11-49f7-83a7-c58800fe130f
# ╟─de63b6f5-b534-4f3e-8bd6-149bf8738ef3
# ╟─216b21f9-38f3-4ce0-a648-8d9a1171b4b6
# ╟─5aad199f-be14-4b9d-b1f0-1e13f553f6fc
