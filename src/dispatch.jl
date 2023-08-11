struct StickySymTuple{X,N}
    data::NamedTuple{X,NTuple{N,Sym}}
end
Base.getproperty(mnt::StickySymTuple, name::Symbol) = getproperty(getfield(mnt, :data), name)
Base.getindex(mnt::StickySymTuple, i::Int) = getfield(mnt, :data)[i]

struct cosHold{T,X}
    angle::T
    expession::X
end

# redefine cosθij
function cosθij(k, σs::StickySymTuple{(:σ1, :σ2, :σ3),3}, msq)
    i, j, _ = ijk(k)
    θ = SymPy.symbols(Symbol("θ_", i, j), real=true)
    cosHold(θ, cosθij(k, getfield(σs, :data), msq))
end


# rederine cosζ
label(wr::WignerRotation) = "^$(wr.k)_" *
                            (ispositive(wr) ? "+" : "-") *
                            (iseven(wr) ? "e" : "0")
#
for N in (0, 2, 3)
    eval(:(
        function cosζ(wr::WignerRotation{$(N)},
            σs::StickySymTuple{(:σ1, :σ2, :σ3),3}, msq)
            ζ = SymPy.symbols("ζ" * label(wr), real=true)
            return cosHold(ζ, cosζ(wr, getfield(σs, :data), msq))
        end
    ))
end

# redefine wignerd_doublearg
function wignerd_doublearg(two_j, two_λ1, two_λ2, cosθ::cosHold)
    half = 1 / Sym(2)
    (abs(two_λ1) > two_j || abs(two_λ2) > two_j) && return zero(cosθ.angle)
    return WignerD(two_j * half, two_λ1 * half, two_λ2 * half,
        0, cosθ.angle, 0)
end


# redefine Clebsches
function CG_doublearg(two_j1, two_m1::Sym, two_j2, two_m2::Sym, two_j, two_m::Sym)
    half = 1 / Sym(2)
    CG(two_j1 * half, two_m1 * half,
        two_j2 * half, two_m2 * half,
        two_j * half, two_m * half)
end