include("./modules/LC_verlet.jl")
include("./modules/basis_naive.jl")
include("./modules/basis_improved.jl")
include("./modules/Particle.jl")
using Plots, LinearAlgebra, .LC_verlet, .basis_naive, .basis_improved

### PARAMS ###
σ, ϵ, L, r_cut, m = 1, 1, 25, 2.5, 1
nc = div(L, r_cut, RoundDown) |> Int
##############

function uniform_sampler(a, b)
    return rand(2,1)*(b-a) .+ a
end

particles = []

for i = 1:10
    push!(particles, Particle.p(m, uniform_sampler(5,20), uniform_sampler(-1,1), [0,0], [0,0]))
end


x, y, e, u = @time LC_verlet.run(particles, 1e-2, 10,r_cut, σ, ϵ, L, nc)

particles

plot(e, legend=:outerleft, label="Kin E")
plot!(u, label="Pot E")
plot!((e+u), label="Total Energy")
savefig("./images/10particles_energy_LC.pdf")

scatter(x, y,title="N=10, t:0->10, dt:1e-2")
savefig("./images/10particles_trajectory_basis.pdf")

