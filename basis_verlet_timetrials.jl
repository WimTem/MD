include("./modules/basis_verlet.jl")
include("./modules/basis_verlet2.jl")
include("./modules/LC_verlet.jl")
include("./modules/Particle.jl")
using .basis_verlet, .basis_verlet2, .LC_verlet, .Particle, Plots, LinearAlgebra

### PARAMS ###
σ, ϵ, L, r_cut = 1, 1, 25, 2.5
nc = Int64(div(L, r_cut, RoundDown))
##############

##Sample from (-1, 1)
function uniform_sampler(a, b)
    return rand(2,1)*(b-a) .+ a
end

##Generate 10 random particles, medium speed.
set = Array{Particle.p}(undef, 50)
for i = 1:50
    set[i] = Particle.p(1,2*uniform_sampler(-1,1).+[13,13], uniform_sampler(-1, 1), [0,0], [0,0])
end

set2 = Array{Particle.p}(undef, 50)
for i = 1:50
    set2[i] = Particle.p(1,2*uniform_sampler(-1,1).+[13,13], uniform_sampler(-1, 1), [0,0], [0,0])
end

#O(n^2)
x_basis, y_basis = @time basis_verlet.run(set, 1e-3, 1, r_cut, σ, ϵ, L)

#O(n^2/2)
x_basis2, y_basis2 = @time basis_verlet2.run(set2, 1e-3, 1, r_cut, σ, ϵ, L)
