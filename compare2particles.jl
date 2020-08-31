include("./modules/LC_verlet.jl")
include("./modules/basis_naive.jl")
include("./modules/basis_improved.jl")
include("./modules/Particle.jl")
using Plots, LinearAlgebra, .LC_verlet, .basis_naive, .basis_improved

### PARAMS ###
σ, ϵ, L, r_cut, m = 1, 1, 25, 2.5, 1
nc = div(L, r_cut, RoundDown) |> Int
##############

p1 = Particle.p(m, [3,3], [2,4], [0,0], [0,0])
p2 = Particle.p(m, [3,4+2^(1/6)], [2,3], [0,0], [0,0])

x1, y1, e1 = @time LC_verlet.run([p1, p2], 1e-3, 10, r_cut, σ, ϵ, L, nc)

p1 = Particle.p(m, [3,3], [2,4], [0,0], [0,0])
p2 = Particle.p(m, [3,4+2^(1/6)], [2,3], [0,0], [0,0])

x2, y2, e2 = @time basis_improved.run([p1, p2], 1e-3, 10, r_cut, σ, ϵ, L)

scatter(x1[1:100:end],y1[1:100:end], label="LC", xlims=(0,25), ylims=(0,25))
scatter!(x2[1:100:end],y2[1:100:end], label="Basis_Improved")
savefig("./images/2particles_trajectory.pdf")

plot(e1[1:1000:end], label="LC", ylims=(0,20), legend=:outerleft, title="t: 0->10, dt: 1e-3", xaxis="t[s]", yaxis="E[J]")
plot!(e2[1:1000:end], label="Basis")
savefig("./images/2particles_energy.pdf")