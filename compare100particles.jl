include("./modules/basis_improved.jl")
include("./modules/Particle.jl")
using Plots, LinearAlgebra, .basis_improved

### PARAMS ###
σ, ϵ, L, r_cut, m = 1, 100, 25, 2.5, 1
nc = div(L, r_cut, RoundDown) |> Int
##############

function uniform_sampler(a, b)
    return rand(2,1)*(b-a) .+ a
end

uniform_sampler(1,2)


particles = []
for i = 1:22
    push!(particles, Particle.p(m, [i*2^(1/6) , 1], [0,0], [0,0],[0,0]))
    push!(particles, Particle.p(m, [i*2^(1/6) , 1+2^(1/6)], [0,0], [0,0],[0,0]))
    push!(particles, Particle.p(m, [i*2^(1/6) , 1+2*2^(1/6)], [0,0], [0,0],[0,0]))
end

push!(particles, Particle.p(m, [1,10], [3, -3], [0,0], [0,0]))


x, y, e, u = @time basis_improved.run(particles, 1e-2, 10, r_cut, σ, ϵ, L)

anim = @animate for i ∈ 1:10:1e3
    scatter(x[Int(i),:], y[Int(i),:], xlims=(0,25), ylims=(0,25))
end
gif(anim, "./images/anim_fps15_100particles.gif", fps = 15)


plot(e[1:100:end], legend=:outerleft, label="Kin E")
plot!(u[1:100:end], label="Pot E")
plot!(e[1:100:end]+u[1:100:end], label="Total Energy")


