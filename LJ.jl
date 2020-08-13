using Plots, LinearAlgebra

####################################
## Enter parameters
ϵ, σ = 1, 1
r_cut = 2.5σ
####################################

r = .5:0.01:8
function E_rep(x)
    return 4*ϵ*(σ/x)^12
end

function E_attr(x)
    return -4*ϵ*(σ/x)^6
end

function LJ(x)
    return 4*ϵ*((σ/x)^12 - (σ/x)^6)
end

plot(r, LJ.(r), xlims=(0,3), ylims=(-2,2), lw=2, label="Lennard-Jones Potential", xlabel="r[m]", ylabel="V[J]")
plot!(r, E_rep.(r), linestyle=:dash, label="Repulsive Term")
plot!(r, E_attr.(r), linestyle=:dash, label="Attractive Term")
plot!(0:.1:1.5, [-1 for i in 0:.1:1.5], linestyle=:dot, label="ϵ")
plot!([2^(1/6) for i in 0:-.1:-1], [i for i in 0:-.1:-1], linestyle=:dot, label="Equilibrium: 2^(1/6)σ")
scatter!([1], [0], label="σ", markersize=3)

savefig("images/LJ_pot.pdf")