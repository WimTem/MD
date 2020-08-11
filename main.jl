using Plots, LinearAlgebra, Random

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
ϵ, σ = 1, 1
L, r_cut = 25, 2.5σ
N = 2
t_end = 10
dt = 1e-2
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


x = init_x(N)
v = init_v(N)


#Initialize Grid
#append!(result[1,1], [Particle(m, [0,0], x[:,i], v[:,i], [0,0]) for i in 1:N])

append!(result[5,1], [Particle(m, [1, 11], [5, 0], [0, 0], [0,0])])

append!(result[5,10], [Particle(m, [24, 11], [-5, 0], [0, 0], [0,0])])

result

# Calculate new F in neighourhood 
function compF_LC(grid, nc, r_cut)
    for i = 1:nc
        for j = 1:nc
            for k in grid[i, j]
                k.F = [0, 0]
            end
        end
    end
    for i = 1:nc
        for j = 1:nc
            for k in grid[i, j]
                for l = i-1:i+1
                    for m = j-1:j+1
                        if m > 0 && m <= nc
                            if l > 0 && l <= nc
                                for n in grid[l, m]
                                    if k != n
                                        force(k, n, r_cut)
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
            for m in grid[i, j]
                for k = 1:2
                    if m.x[k] < 0 || m.x[k] > L
                        m.v[k] = -m.x[k]
                    end
                end
            end
        end
    end
end


function force(p1, p2, r_cut)
    r_vector = p2.x - p1.x
    rmag = norm(p1.x - p2.x)
    if rmag <= r_cut
        s = (σ/rmag)^6
        f = 24ϵ*(s/rmag)*(1-2s)*r_vector
        p1.F = p1.F + f
    end
end


# N = amount of particles
# n = amount of datapoints
function timeIntegration_LC(t::Real, dt::Real, t_end::Real, grid, n::Int64)
    result_x = zeros(N,n+1)
    result_y = zeros(N,n+1)
    compF_LC(grid, nc, r_cut)
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
                if length(grid[i, j]) > 0
                    for k in grid[i, j]
                        z += 1
                        result_x[z, count] = k.x[1]
                        result_y[z, count] = k.x[2]
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

plot(x[1,200:240], y[1,200:240], lw=3, linestyle=:dashdot)
plot!(x[2,200:240], y[2,200:240], lw=3, linestyle=:dashdot)

