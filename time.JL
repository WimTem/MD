include("./modules/basis_verlet.jl")
include("./modules/LC_verlet.jl")
include("./modules/Particle.jl")
using .basis_verlet, .LC_verlet, .Particle, Plots, LinearAlgebra

p1 = Particle.Part(1, [1, 1], [0, 0], [0, 0], [0,0])
p2 = Particle.Part(1, [1, 2], [0, 0], [0, 0], [0,0])

### PARAMS ###
σ, ϵ, L, r_cut = 1, 1, 25, 2.5
nc = Int64(div(L, r_cut, RoundDown))
##############

result = Array{Any}(undef, nc,nc)
    for i = 1:nc
        for j = 1:nc
            result[i, j] = Array{Particle.Part}(undef, 0)
        end
    end
result
    append!(result[1,1], [p1])
    result
    append!(result[1,1], [p2])


result

x_LC, y_LC = @time LC_verlet.run(result, .1, 10, 2.5, 1, 1, L, nc)

p1 = Particle.Part(1, [1, 1], [0, 0], [0, 0], [0,0])
p2 = Particle.Part(1, [1, 2], [0, 0], [0, 0], [0,0])
x_basis, y_basis = @time basis_verlet.run([p1, p2], 0.1, 10, 2.5, 1, 1, 25)



scatter(y_basis, label="Base method", layout=(2,1), legend=:outertopleft)
scatter!(y_LC, label="LC method", layout=(2,1), legend=:outertopleft)

savefig("Indrukwekkende nauwkeurigheid.pdf")