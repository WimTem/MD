using Plots, LinearAlgebra, Random, Animations

mutable struct Particle
    m
    x
    v
    F
    F_old
    Particle(m, x, v, F, F_old) = new(m, x, v, F, F_old)
end

# Divide total domain Ω in subdomain cells, length of square cell = l
# l >= r_cut : Interactions of particles limited to cell + neighbouring cells
# ∴ cell = rcut x rcut
# Fᵢ ≈ ∑(kc∈ N(ic)) ∑j ∈ (particles in kc) Fᵢⱼ

# Assume square domain (0,L)x(0,L)
# Max number of cells per dim: nc = RoundDown(L/r_cut)
# N: number of particles

####################################
## Enter parameters
m = 1
ϵ, σ = 10, 1
L, r_cut = 25, 2.5σ
N = 4
t_end = 100
dt = 1e-3
l = r_cut
####################################
nc = Int64(div(L, r_cut, RoundDown))

function uniform_sampler(a, b)
    return rand(2,1)*(b-a) .+ a
end

function init_x(n)
    result = zeros(2,n)
    for i = 1:n
        result[:,i] = uniform_sampler(0, 2.5)
    end
    return result
end

function init_v(n)
    result = zeros(2,n)
    for i = 1:n
        result[:, i] = uniform_sampler(-5, 5)
    end
    return result
end

result = Array{Any}(undef, nc,nc)
for i = 1:nc
    for j = 1:nc
        result[i, j] = Array{Particle}(undef, 0)
    end
end

#Initialize Grid
#append!(result[1,1], [Particle(m, [0,0], x[:,i], v[:,i], [0,0], [0, 0]) for i in 1:N])

append!(result[1, 1], [Particle(m, [12,12], [0, 0], [0, 0], [0,0])])
append!(result[1,1], [Particle(m, [12, 13], [0, 0], [0, 0], [0,0])])
append!(result[1,1], [Particle(m, [14, 2], [-10, 0], [0, 0], [0,0])])
append!(result[1,1], [Particle(m, [4,2], [10, 0], [0, 0], [0,0])])

result

# Calculate new F in neighourhood 
function compF_LC(grid, nc, r_cut)
    #Loop over cells
    for i = 1:nc
        for j = 1:nc
            #Loop over elements in cell
            if length(grid[i, j]) > 0
                #Loop over neighbourhood cells of each element
                for el1 in grid[i, j]
                    for k in i-1:i+1
                        for l in j-1:j+1
                            if k > 0 && k <= nc
                                if l > 0 && l <= nc
                                    for el2 in grid[k, l]
                                        if el1 != el2
                                            force(el1, el2, r_cut)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function compX_LC(grid, nc, l, dt)
    for i = 1:nc
        for j = 1:nc
            for k in grid[i, j]
                updateX(k, dt)
            end
        end
    end
    moveParticles_LC(grid, nc, l)
end

function compV_LC(grid, nc, l ,dt)
    for i = 1:nc
        for j = 1:nc
            for k in grid[i, j]
                updateV(k, dt)
            end
        end
    end
end

function updateX(p, dt)
    a = dt*0.5/p.m
    p.x = p.x + dt*(p.v + a*p.F)
    p.F_old = p.F
end

function updateV(p, dt)
    a = dt*0.5/p.m
    p.v = p.v + a*(p.F + p.F_old)
end

function moveParticles_LC(grid, nc, l)
    for i = 1:nc
        for j = 1:nc
            ## Reflecting boundaries
            for k in grid[i, j]
                for m = 1:2
                    if k.x[m] <= 0
                        k.v[m] = -k.v[m]
                        k.x[m] = 1e-3 #Small positive required, see lline 147
                    end
                    if k.x[m] > L
                        k.v[m] = -k.v[m]
                        k.x[m] = L
                    end
                end
            end
            ## Move back to correct cell
            for k in grid[i, j]
                true_col = Int(div(k.x[1], r_cut, RoundUp))
                true_row = Int(div(k.x[2], r_cut, RoundUp))
                if true_col != j || true_row != i
                    push!(grid[true_row, true_col], k)
                    filter!(e->e!=k, grid[i, j])
                end
            end
        end
    end
end

function force(p1, p2, r_cut)
    rmag = norm(p1.x - p2.x)
    r_vector = (p2.x - p1.x)/rmag ## Normalized vector
    if rmag <= r_cut
        s = (σ/rmag)^6
        p1.F = 24ϵ*(s)*(1-2s)*r_vector
    end
end


# N = amount of particles
# n = amount of datapoints
function timeIntegration_LC(t::Real, dt::Real, t_end::Real, grid, n::Int64)
    result_x = zeros(n, N)
    result_y = zeros(n, N)
    count = 0
    while (t < t_end)
        count += 1
        t += dt
        z = 0
        compX_LC(grid, nc, l, dt)
        compF_LC(grid, nc, r_cut)
        compV_LC(grid, nc, l ,dt)
        for i = 1:nc
            for j = 1:nc
                if length(grid[i, j]) > 0 && count < n
                    for k in grid[i, j]
                        z += 1
                        result_x[count, z] = k.x[1]
                        result_y[count, z] = k.x[2]
                    end
                end
            end
        end
    end
    return result_x, result_y
end

function main(grid, dt, t_end)
    n = Int(t_end/dt)
    result_x, result_y = timeIntegration_LC(0, dt, t_end, grid, n)
    return result_x, result_y
end

x, y = main(result, dt, t_end)

x

scatter([x[1,1]], [y[1,1]],  xlims=(0,25), ylims=(0,25), legend=false, color=:red, markersize=5)
scatter!([x[1,2]], [y[1,2]],  xlims=(0,25), ylims=(0,25), legend=false, color=:blue, markersize=5)


anim = @animate for i = 1:10000
    scatter(x[i,:], y[i,:],  xlims=(0,25), ylims=(0,25), legend=false, palette=:blues, markersize=5)
end every 1000

gif(anim, "images/evenwicht_realtime.gif", fps=24)

