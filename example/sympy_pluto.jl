### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 89fb4510-36f4-11ee-09a5-b95a60001661
begin
	cd(joinpath(@__DIR__))
	using Pkg
	Pkg.activate("..")
	
	using ThreeBodyDecay
	using SymPy
	import SymPy.PyCall
	# 
	PyCall.pyimport_conda("sympy.physics.wigner", "sympy")
	import_from(sympy.physics.wigner)
	# wigner_d_small(Sym(1),θ)[1,1] # = cos(θ) d^1_00(θ)
	# 
	PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
	import_from(sympy.physics.quantum.spin, (:WignerD,), typ=:Any)
	# WignerD
end

# ╔═╡ 37eca537-c92f-4006-8f45-a15ca3174f9e
md"""
# ThreeBodyDecay with SymPy

The notebook calls `ThreeBodyDecay` implementation passing `SymPy` object.
With redefinion of a few function, the amplitude function spits nice symbolic expredssion. Nice!
"""

# ╔═╡ f8f7d1d6-a563-4226-a359-d9f5da933fc8
begin
	@syms θ z
	# 
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

# ╔═╡ 35a8d6f3-81a7-4036-ac6f-7d42975db0bf
begin # Ξc(J) -> Ξc(JP) [-> Ξc(1/2+) π(0+)] π(0+) 
	reaction = "3/2-" => ("1/2+", "0-", "0-")
	js, Ps = reaction |> spinparity
	tbs = ThreeBodySystem(ms, js)
end

# ╔═╡ c74432fa-a8e8-490e-8b75-5367a2d285f1
begin # resonance
	Rjp = jp"3/2+"
	R(σ) = Sym("R")
end

# ╔═╡ c85ddfd8-0948-4d84-a9f6-caecda46afaa
dc = let
	dcv = DecayChainsLS(3, R; two_s=Rjp.j |> x2, parity=Rjp.p, Ps, tbs)
	dc = dcv[1, 1]
	@show dc.HRk.two_ls
	@show dc.Hij.two_ls
	dc
end

# ╔═╡ 228c0817-a067-4ce6-b958-34a636228011
md"""
## Functionality extension via magical dispatch
The key is to redefine `wignerd_doublearg` to give `SymPy.WignerD` instead of numerical `PartialWaveFunctions.wignerd`
"""

# ╔═╡ a0c25e13-3008-45ff-ae73-858edab939d4
begin
	import ThreeBodyDecay: cosθij, cosζ
	import ThreeBodyDecay.PartialWaveFunctions: wignerd_doublearg
	import ThreeBodyDecay: MandestamTuple, WignerRotation, wr
	
	struct cosHold{T}
	    angle::T
	end	
	function cosθij(k, σs::MandestamTuple{Sym}, msq)
	    i, j, _ = ijk(k)
	    θ = Sym(Symbol("θ_", i, j)) |> first
	    cosHold(θ)
	end
	# 
	label(wr::WignerRotation) = "^$(wr.k)_" *
	                            (ispositive(wr) ? "+" : "-") *
	                            (iseven(wr) ? "e" : "0")
	# 
	cosζ(wr::WignerRotation{0}, σs::MandestamTuple{Sym}, msq) = Sym("ζ" * label(wr)) |> cosHold
	cosζ(wr::WignerRotation{2}, σs::MandestamTuple{Sym}, msq) = Sym("ζ" * label(wr)) |> cosHold
	cosζ(wr::WignerRotation{3}, σs::MandestamTuple{Sym}, msq) = Sym("ζ" * label(wr)) |> cosHold
	# 
	function wignerd_doublearg(two_j, two_λ1, two_λ2, cosθ::cosHold)
	    half = 1 / Sym(2)
	    WignerD(two_j * half, two_λ1 * half, two_λ2 * half,
	        0, cosθ.angle, 0)
	end
end

# ╔═╡ 8e7a77fa-ed71-4295-92dd-a65df648b3f7
[amplitude(dc, σs, two_λs; refζs=(3, 1, 1, 3)).doit() |> simplify 
	for two_λs in itr(js)] |> vec

# ╔═╡ 52eb2a1b-64f5-4240-a6e8-bb9d99afc15b
sum(itr(js)) do two_λs
	amplitude(dc, σs, two_λs; refζs=(3, 1, 1, 3)).doit()^2
end |> simplify 

# ╔═╡ Cell order:
# ╟─37eca537-c92f-4006-8f45-a15ca3174f9e
# ╠═89fb4510-36f4-11ee-09a5-b95a60001661
# ╠═f8f7d1d6-a563-4226-a359-d9f5da933fc8
# ╠═ebcd010a-ea59-4ef7-8d49-99c263ce20ac
# ╠═cb54532e-c0d5-46d1-9eb9-30b63a72f5d5
# ╠═864208ca-42a5-40ff-a0c2-f7fccc4f6cbe
# ╠═35a8d6f3-81a7-4036-ac6f-7d42975db0bf
# ╠═c74432fa-a8e8-490e-8b75-5367a2d285f1
# ╠═c85ddfd8-0948-4d84-a9f6-caecda46afaa
# ╟─228c0817-a067-4ce6-b958-34a636228011
# ╠═a0c25e13-3008-45ff-ae73-858edab939d4
# ╠═8e7a77fa-ed71-4295-92dd-a65df648b3f7
# ╠═52eb2a1b-64f5-4240-a6e8-bb9d99afc15b
