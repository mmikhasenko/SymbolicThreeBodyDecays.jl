using SymbolicThreeBodyDecays
using SymbolicThreeBodyDecays.ThreeBodyDecay
using SymbolicThreeBodyDecays.SymPy
using Test

@syms m1 m2 m3 m0
@syms σ1 σ2 σ3
# 
ms = ThreeBodyMasses(m1, m2, m3; m0)
σs = Invariants(ms; σ1, σ2)


@testset "Properties of StickySymTuple" begin
    sσs = StickySymTuple(σs)
    @test sσs.σ1 == σs.σ1 && sσs.σ2 == σs.σ2 && sσs.σ3 == σs.σ3
    @test sσs[1] == σs[1] && sσs[2] == σs[2] && sσs[3] == σs[3]
end

@testset "Evaluation with StickySymTuple gives cosHold" begin
    @test cosθ12(σs, ms^2) isa Sym
    @test cosθ12(σs |> StickySymTuple, ms^2) isa cosHold

    @test cosζ(wr(1, 2, 1), σs, ms^2) isa Sym
    @test cosζ(wr(1, 2, 1), σs |> StickySymTuple, ms^2) isa cosHold
end


function spinparity(p)
    pt = (p[2]..., p[1])
    jpv = str2jp.(pt)
    getproperty.(jpv, :j) .|> x2 |> ThreeBodyDecay.SpinTuple,
    getproperty.(jpv, :j) .|> x2 |> ThreeBodyDecay.ParityTuple
end

dc = let
    # Ξc(J) -> Ξc(3/2+) [-> Ξc(1/2+) π(0+)] π(0+) 
    reaction = "1/2+" => ("1/2+", "0-", "0-")
    js, Ps = reaction |> spinparity
    tbs = ThreeBodySystem(ms, js)
    #
    Rjp = jp"3/2+"
    # 
    DecayChainLS(3, σ -> Sym("R"); two_s=Rjp.j |> x2, parity=Rjp.p, Ps, tbs)
end

expr = amplitude(dc, σs |> StickySymTuple, (1, 0, 0, 1); refζs=(1, 2, 3, 1)).doit() |> simplify
expr_symbols = expr.free_symbols |> collect

@testset "Simplification with StickyTuple" begin
    @test σ1 ∉ expr_symbols
    @test σ2 ∉ expr_symbols
    @test σ3 ∉ expr_symbols
    @test sympy.symbols("θ_12", real=true) ∈ expr_symbols
end

@testset "Labels for Wigner Rotations" begin
    triv = SymbolicThreeBodyDecays.label(TriavialWignerRotation(1))
    for (i,j,k) in Iterators.product(1:3,1:3,1:3)
        p = SymbolicThreeBodyDecays.label(wr(i,j,k)) => 
            ((i != j) ? "^$(k)_$(i)($(j))" : triv)
        @test p[1] == p[2]
    end
end


@testset "Free zeta angles in amplitude" begin
    free_symbols(refζs) = amplitude(dc, σs |> StickySymTuple, (1, 0, 0, 1);
        refζs).doit().free_symbols .|> string |> Set
    # 
    @test free_symbols((3,1,1,3)) == Set(["θ_12", "R"])
    @test free_symbols((2,1,1,3)) == Set(["θ_12", "R", "ζ^1_3(2)"])
    @test free_symbols((3,1,1,1)) == Set(["θ_12", "R", "ζ^0_3(1)"])
    @test free_symbols((1,1,1,2)) == Set(["θ_12", "R", "ζ^0_3(2)", "ζ^1_3(1)"])
end
    
