### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 89fb4510-36f4-11ee-09a5-b95a60001661
# ╠═╡ show_logs = false
begin
    cd(joinpath(@__DIR__, ".."))
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()
    Pkg.add("Plots")
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
# Lb decays with open flavor (Pc search)

The notebook computes the angular distributions in the decay

$\varLambda_b^0 \to D_{s1}^- (\to D^{*0} K^-)\, \varLambda_c^+$

There are two possible LS couplings in the decays of $D_{s1}$, the S and D waves,
and several possible LS in the $\Lambda_b^0$ decays due to unconstrained parity in the weak-current process.
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
# Lb(1/2) -> Ds*(1+) [-> D*(1-) K(0-)] Lc(1/2+) 
# 
function decay_via_xic2645(Rjp; P0='+')
    # 
    ifstate = "1/2"*P0 => ("1-", "0-", "1/2+")
    js, Ps = ifstate |> spinparity
    tbs = ThreeBodySystem(ms, js)
    # 
    # resonance
    (; j, p) = str2jp(Rjp)
    R(σ) = Sym("R")
    # 
    # decay chain
    dc_all = DecayChainsLS(3, R; two_s=j |> x2, parity=p, Ps, tbs)
    return dc_all
end

# ╔═╡ 481b0408-cc1d-4d9c-b3bf-a86f5eb21b0d
function I(dc)
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
I.(decay_via_xic2645("1+", P0='+')) |> vec

# ╔═╡ 0851bca3-6d67-4070-a449-5a6450b18078
I.(decay_via_xic2645("1+", P0='-')) |> vec

# ╔═╡ 2a01b8da-3066-469e-a67d-3738816f6c79
begin
	plot(title="\$\\Lambda_b^0 \\to D_{s1}^- (\\to D^{*0} K^-)\\, \\Lambda_c^+\$")
	for P0 in ['+','-']
		dcv = decay_via_xic2645("1+"; P0) |> vec
		@syms θ12::real=>"θ_12" cosθ
		map(dcv) do dc
			I_jp = I(dc)
			# 
			L, S = div.(dc.HRk.two_ls,2)
			l, s = div.(dc.Hij.two_ls,2)
			# 
	        expr = I_jp.subs(Dict(
				dc.Xlineshape(0) => Sym(1),
				θ12 => acos(cosθ)
			), simultaneous=true) + 1e-7 * cosθ + 0.003*(L+S+l)
			plot!(expr, -1, 1; ls = (P0 == '+' ? :solid : :dash),
				lab=(P0 == '+' ? "PC" : "PV") * ", $(L) $(S) $(l)", lw=1.5, c=(S+l))
		end        
	end
	plot!(ylim=(0,:auto), legend_title="\$\\Lambda_b^0\$, L S l",
			xlab="cosine of \$D_{s1}^*\$ helicity angle")
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
# ╠═481b0408-cc1d-4d9c-b3bf-a86f5eb21b0d
# ╟─9c2366cd-1a54-4d6f-b00f-591871c9d794
# ╠═6c2d9cc6-38f9-47c2-aa83-0986394ed58b
# ╠═0851bca3-6d67-4070-a449-5a6450b18078
# ╠═2a01b8da-3066-469e-a67d-3738816f6c79
