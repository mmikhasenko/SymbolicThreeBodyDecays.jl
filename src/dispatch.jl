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
    θ = Sym(Symbol("θ_", i, j))[1]
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
            ζ = Sym("ζ" * label(wr))
            return cosHold(ζ, cosζ(wr, getfield(σs, :data), msq))
        end
    ))
end

# redefine wignerd_doublearg
function wignerd_doublearg(two_j, two_λ1, two_λ2, cosθ::cosHold)
    half = 1 / Sym(2)
    WignerD(two_j * half, two_λ1 * half, two_λ2 * half,
        0, cosθ.angle, 0)
end


