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
    import Plots.PlotMeasures: mm
end

# ╔═╡ 37eca537-c92f-4006-8f45-a15ca3174f9e
md"""
# Studies of Omegac** spin

The calculations revieals the angular distribution from the [LHCb paper](https://inspirehep.net/literature/1879440).
It is a weak decay of $\Omega_c^0$ baryon

$\varOmega_c^-(1/2^\pm) \to \varOmega_c^{**0}(J^P) \pi^-$

The $\varOmega_c^{**0}(J^P)$ resonances appear in the $\varLambda_c^+ K^-$ system.
The angular distributions are not sensitive to the parity. Still a non-trivial angular dependence appears since $\Omega_b^-$ spin only avarages out $\pm 1/2$ spin projections. 
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
# Ωb (1/2+) -> Ωc(JP) [-> Lc(1/2+) K(0+)] π(0+) 
# 
function decay_via_xic2645(jpR, P0='+')
    # 
    ifstate = "1/2"*P0 => ("1/2+", "0-", "0-")
    js, Ps = ifstate |> spinparity
    tbs = ThreeBodySystem(ms, js)
    # 
    # resonance
    Rjp = str2jp(jpR)
    R(σ) = Sym("R")
    # 
    # decay chain
    dc_all = DecayChainsLS(3, R; two_s=Rjp.j |> x2, parity=Rjp.p, Ps, tbs)
	size(dc_all) != (1,1) && @info "More couplings"
    return dc_all |> first
end

# ╔═╡ a0c25e13-3008-45ff-ae73-858edab939d4
dc_test = let
    dc = decay_via_xic2645("1/2+")
    @show dc.HRk.two_ls
    @show dc.Hij.two_ls
    dc
end;

# ╔═╡ 481b0408-cc1d-4d9c-b3bf-a86f5eb21b0d
function I(jp0, P0='+')
    dc = decay_via_xic2645(jp0, P0)
    full_amplitude = sum(itr(dc.tbs.two_js)) do two_λs
        A = amplitude(dc, σs |> StickySymTuple, two_λs .|> Sym; refζs=(3, 1, 1, 3)).doit()
		abs2(A)
    end
    full_amplitude |> simplify |> simplify
end

# ╔═╡ 9c2366cd-1a54-4d6f-b00f-591871c9d794
md"""
## Test different $J^P$ of the resonance
"""

# ╔═╡ 6c2d9cc6-38f9-47c2-aa83-0986394ed58b
I("1/2+",'+') == I("1/2-",'+') == I("1/2+",'-') == I("1/2-",'-') && I("1/2+")

# ╔═╡ 0edd1ecd-3a77-4b54-bb27-de20207e3ca1
I("3/2+",'+') == I("3/2-",'+') == I("3/2+",'-') == I("3/2-",'-') && I("3/2+")

# ╔═╡ 7ed77d8f-eee6-4ccf-ae14-005d3529adab
I("5/2+",'+') == I("5/2-",'+') == I("5/2+",'-') == I("5/2-",'-') && I("5/2+")

# ╔═╡ 2a01b8da-3066-469e-a67d-3738816f6c79
begin
    plot(title="\$\\varOmega_b^-(1/2^\\pm) \\to \\varOmega_c^{**}(J^P) \\pi^-\$")
    for jp in vec(string.(1:2:5) .* "/2+")
		I_jp = I(jp)
		@syms θ12::real=>"θ_12" cosθ
        expr = I_jp.subs(Dict(
			dc_test.Xlineshape(0) => Sym(1),
			θ12 => acos(cosθ)
		), simultaneous=true) + 1e-7 * cosθ
        plot!(expr, -1, 1, lab=jp, lw=1.5)
    end
    plot!(ylim=(0, :auto), legend_title="JP",
		xlab="cosine of \$\\varOmega_c^{**0}\$ helicity angle")
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
# ╠═0edd1ecd-3a77-4b54-bb27-de20207e3ca1
# ╠═7ed77d8f-eee6-4ccf-ae14-005d3529adab
# ╠═2a01b8da-3066-469e-a67d-3738816f6c79
