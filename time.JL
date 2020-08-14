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


scatter(x_LC, y_LC, legend=:outertopleft, label="LC methode")
scatter!(x_basis, y_basis, legend=:outertopleft, label="Basis methode")
savefig("trajectories.pdf")

for i in []
    print("ja")
end

savefig("topresultaat.pdf")

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

