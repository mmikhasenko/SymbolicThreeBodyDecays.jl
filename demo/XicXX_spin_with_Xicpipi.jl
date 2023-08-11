### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 89fb4510-36f4-11ee-09a5-b95a60001661
# ╠═╡ show_logs = false
begin
    cd(joinpath(@__DIR__, ".."))
    using Pkg
    Pkg.activate(".")
    # 
    using SymbolicThreeBodyDecays
    using ThreeBodyDecay
    using SymPy
	# 
	using Plots
	import Plots.PlotMeasures:mm
end

# ╔═╡ 37eca537-c92f-4006-8f45-a15ca3174f9e
md"""
# ThreeBodyDecay with SymPy

The notebook calls `ThreeBodyDecay` implementation passing `SymPy` object.
With redefinion of a few function, the amplitude function spits nice symbolic expredssion. Nice!
"""

# ╔═╡ 8b28c3c6-9e94-419d-9e35-ac82edc3a832
theme(:wong2, frame=:box, grid=false, minorticks=true,
	guidefontvalign=:top, guidefonthalign=:right,
	xlim=(:auto, :auto), ylim=(:auto, :auto),
	lw=1.2, lab="", colorbar=false,
	bottom_margin=3mm, left_margin=3mm)

# ╔═╡ f8f7d1d6-a563-4226-a359-d9f5da933fc8
begin
    @syms m1 m2 m3 m0
    @syms σ1 σ2 σ3
end;

# ╔═╡ ebcd010a-ea59-4ef7-8d49-99c263ce20ac
ms = ThreeBodyMasses(m1, m2, m3; m0)

# ╔═╡ cb54532e-c0d5-46d1-9eb9-30b63a72f5d5
σs = Invariants(ms; σ1, σ2)

# ╔═╡ 864208ca-42a5-40ff-a0c2-f7fccc4f6cbe
function spinparity(p)
    pt = (p[2]..., p[1])
    jpv = str2jp.(pt)
    getproperty.(jpv, :j) .|> x2 |> ThreeBodyDecay.SpinTuple,
    getproperty.(jpv, :p) |> ThreeBodyDecay.ParityTuple
end

# ╔═╡ c85ddfd8-0948-4d84-a9f6-caecda46afaa
# Ξc(J) -> Ξc(3/2+) [-> Ξc(1/2+) π(0+)] π(0+) 
# 
function decay_via_xic2645(jp0)
    # 
    ifstate = jp0 => ("1/2+", "0-", "0-")
    js, Ps = ifstate |> spinparity
    tbs = ThreeBodySystem(ms, js)
    # 
    # resonance
    Rjp = jp"3/2+"
    R(σ) = Sym("R")
    # 
    # decay chain
    dc = DecayChainLS(3, R; two_s=Rjp.j |> x2, parity=Rjp.p, Ps, tbs)
    return dc
end

# ╔═╡ a0c25e13-3008-45ff-ae73-858edab939d4
dc_12p = let
    dc = decay_via_xic2645("1/2+")
    @show dc.HRk.two_ls
    @show dc.Hij.two_ls
    dc
end;

# ╔═╡ 481b0408-cc1d-4d9c-b3bf-a86f5eb21b0d
function I(jp0)
    dc = decay_via_xic2645(jp0)
    full_amplitude = sum(itr(dc.tbs.two_js)) do two_λs
		A = amplitude(dc, StickySymTuple(σs), Sym.(two_λs); refζs=(3, 1, 1, 3)).doit()
    	abs2(A)
    end
    full_amplitude |> simplify |> simplify
end

# ╔═╡ 9c2366cd-1a54-4d6f-b00f-591871c9d794
md"""
## Test different $J^P$ of the decay particle
"""

# ╔═╡ 6c2d9cc6-38f9-47c2-aa83-0986394ed58b
I("1/2+")

# ╔═╡ 7763a2ec-a292-4b11-886d-0249b9f45798
I("1/2-")

# ╔═╡ 0edd1ecd-3a77-4b54-bb27-de20207e3ca1
I("3/2+")

# ╔═╡ 7ed77d8f-eee6-4ccf-ae14-005d3529adab
I("3/2-")

# ╔═╡ ebb97af7-0968-4f1f-a63c-7bb1e37eec7a
I("5/2+")

# ╔═╡ 152cd4bb-7a18-4c06-a284-62b158786cf9
I("5/2-")

# ╔═╡ 2a01b8da-3066-469e-a67d-3738816f6c79
begin
	plot(title="Ξc(JP) → Ξc(3/2+)π")
	for jp in vec(string.(1:2:5) .* "/2" .* ['+' '-'])
		expr = I(jp).subs(dc_12p.Xlineshape(0), Sym(1), simultaneous=true)
		@syms θ12::real=>"θ_12" cosθ
		e = expr.subs(θ12, acos(cosθ)) + 1e-7*cosθ
		# @show e
		plot!(e, -1, 1, lab=jp, lw=1.5)
	end
	plot!(ylim=(0,:auto), legend_title="JP")
end

# ╔═╡ Cell order:
# ╟─37eca537-c92f-4006-8f45-a15ca3174f9e
# ╠═89fb4510-36f4-11ee-09a5-b95a60001661
# ╠═8b28c3c6-9e94-419d-9e35-ac82edc3a832
# ╠═f8f7d1d6-a563-4226-a359-d9f5da933fc8
# ╠═ebcd010a-ea59-4ef7-8d49-99c263ce20ac
# ╠═cb54532e-c0d5-46d1-9eb9-30b63a72f5d5
# ╠═864208ca-42a5-40ff-a0c2-f7fccc4f6cbe
# ╠═c85ddfd8-0948-4d84-a9f6-caecda46afaa
# ╠═a0c25e13-3008-45ff-ae73-858edab939d4
# ╠═481b0408-cc1d-4d9c-b3bf-a86f5eb21b0d
# ╟─9c2366cd-1a54-4d6f-b00f-591871c9d794
# ╠═6c2d9cc6-38f9-47c2-aa83-0986394ed58b
# ╠═7763a2ec-a292-4b11-886d-0249b9f45798
# ╠═0edd1ecd-3a77-4b54-bb27-de20207e3ca1
# ╠═7ed77d8f-eee6-4ccf-ae14-005d3529adab
# ╠═ebb97af7-0968-4f1f-a63c-7bb1e37eec7a
# ╠═152cd4bb-7a18-4c06-a284-62b158786cf9
# ╠═2a01b8da-3066-469e-a67d-3738816f6c79
