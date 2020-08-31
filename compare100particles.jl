include("./modules/LC_verlet.jl")
include("./modules/Particle.jl")
using Plots, LinearAlgebra, .LC_verlet

### PARAMS ###
σ, ϵ, L, r_cut, m = 1, 1, 25, 2.5, 1
nc = div(L, r_cut, RoundDown) |> Int
##############

function uniform_sampler(a, b)
    return rand(2,1)*(b-a) .+ a
end

uniform_sampler(1,2)

particles = []

for i = 1:100
    push!(particles, Particle.p(m, uniform_sampler(5,20), uniform_sampler(-5,5), [0,0], [0,0]))
end


x, y, e, u = @time LC_verlet.run(particles, 1, 1e3, r_cut, σ, ϵ, L, nc)

anim = @animate for i ∈ 1:1:1e3
    scatter(x[Int(i),:], y[Int(i),:], xlims=(0,25), ylims=(0,25))
end
gif(anim, "./images/anim_fps15_100particles.gif", fps = 15)

plot(e[1:100:end], legend=:outerleft, label="Kin E")
plot!(u[1:100:end], label="Pot E")
plot!(e[1:100:end]+u[1:100:end], label="Total Energy")


rem(-57,25)