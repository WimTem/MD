include("./modules/basis_verlet2.jl")
include("./modules/Particle.jl")
using Plots, LinearAlgebra, .basis_verlet

### PARAMS ###
σ, ϵ, L, r_cut, m = 1, 1, 25, 2.5, 1
nc = div(L, r_cut, RoundDown) |> Int
##############

p1 = Particle.p(m, [3,3], [2,4], [0,0], [0,0])
p2 = Particle.p(m, [3,4+2^(1/6)], [2,3], [0,0], [0,0])

x, y, e = basis_verlet2.run([p1, p2], 1e-3, 10, r_cut,σ, ϵ, L)

p1 = scatter(x, y, label=["Particle1" "Particle2"], xlims=(0,25), ylims=(0,25), legend=:outertopleft)
p2 = scatter(e[1:100:end], xlabel="t [s]", ylabel="Energy")
scatter(p1, p2, layout=(2,1))
savefig("Fysisch_onmogelijk.pdf")