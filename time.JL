include("./modules/basis_verlet.jl")
include("./modules/LC_verlet.jl")
include("./modules/Particle.jl")
using .basis_verlet, .LC_verlet, .Particle, Plots, LinearAlgebra

### PARAMS ###
σ, ϵ, L, r_cut = 1, 1, 25, 2.5
nc = Int64(div(L, r_cut, RoundDown))
##############

p1 = Particle.Part(1, [1, 1], [1, 10], [0, 0], [0,0])
p2 = Particle.Part(1, [24, 24], [-1, -10], [0, 0], [0,0])

x_LC, y_LC = @time LC_verlet.run([p1, p2], 1e-2, 10, 2.5, 1, 1, L, nc)

p1 = Particle.Part(1, [1, 1], [1, 10], [0, 0], [0,0])
p2 = Particle.Part(1, [24, 24], [-1, -10], [0, 0], [0,0])
x_basis, y_basis = @time basis_verlet.run([p1, p2], 1e-2, 10, 2.5, 1, 1, L)


scatter(x_LC[1:25:end,:], y_LC[1:25:end,:], legend=:outertopleft, label="LC methode")
scatter!(x_basis[10:25:end,:], y_basis[10:25:end,:], legend=:outertopleft, label="Basis methode")
savefig("trajectories.pdf")

p1 = scatter(x_LC[:,1], y_LC[:,1], label="Particle1", xlims=(0,25), ylims=(0,25))
p2 = scatter(x_LC[:,2], y_LC[:,2], label="Particle2", xlims=(0,25), ylims=(0,25))

plot(p1, p2, layout=(2,1))
savefig("uitmuntend.pdf")

result = Array{Any}(undef, nc,nc)
    for i = 1:nc
        for j = 1:nc
            result[i, j] = Array{Particle.Part}(undef, 0)
        end
    end
result
push!(result[1,1], p1)
result[2,1]
append!(result[1,1], [p2])

