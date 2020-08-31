 module LC_verlet
    using LinearAlgebra, Random, Test

    function timeIntegration_LC(t::Real, dt::Real, t_end::Real, grid, N::Int64, n::Int64, r_cut, σ, ϵ, L, nc)
        result_x = zeros(n, N)
        result_y = zeros(n, N)
        result_e = zeros(n)
        count = 0
        for i = 1:n
            count += 1
            compX_LC(grid, dt, r_cut, L, nc)
            compF_LC(grid, σ, ϵ, r_cut, nc)
            compV_LC(grid, dt, nc)
            z = 0
            for i in Iterators.product(1:nc, 1:nc) |> collect
                if length(grid[i[1], i[2]]) > 0 && count <= n
                    for k in grid[i[1], i[2]]
                        z += 1
                        result_x[count, z] = k.x[1]
                        result_y[count, z] = k.x[2]
                        result_e[count] += compoutStatistic_LC(k)
                    end
                end
            end
        end
        return result_x, result_y, result_e
    end

    function compoutStatistic_LC(k)
        v = norm(k.v)
        return .5*k.m*v^2
    end

    function compX_LC(grid, dt, r_cut, L, nc)
        for i in Iterators.product(1:nc, 1:nc) |> collect
            if length(grid[i[1], i[2]]) > 0
                for k in grid[i[1], i[2]]
                    updateX(k, dt)
                end
            end
        end
        moveParticles_LC(grid, nc, L, r_cut)
    end

    function updateX(p, dt)
        p.x = p.x + dt*(p.v + (dt*0.5/p.m)*p.F)
        p.F_old = p.F
    end

    function compF_LC(grid, σ, ϵ, r_cut, nc)
        for i in Iterators.product(1:nc, 1:nc) |> collect
            for p1 in grid[i[1], i[2]]
                for j in neighbour_cells(i[1], i[2])
                    for p2 in grid[j[1], j[2]]
                        if p1 != p2
                            p1.F =  force(p1, p2, σ, ϵ, r_cut)
                        end
                    end
                end
            end
        end
    end

    function moveParticles_LC(grid, nc, L, r_cut)
        for i in Iterators.product(1:nc, 1:nc) |> collect
            for j in grid[i[1], i[2]]
                # Reflecting boundaries
                for m = 1:2
                    if j.x[m] <= 0
                        j.v[m] = -j.v[m]
                        j.x[m] = 1e-6 #Small positive required
                    end
                    if j.x[m] > L 
                        j.v[m] = -j.v[m]
                        j.x[m] = L
                    end
                end
            end
            ## Move back to correct cell
            for j in grid[i[1], i[2]]
                index = div.(j.x, r_cut, RoundUp) .|> Int
                if i[1] != index[1]
                    push!(grid[index[1], i[2]], splice!(grid[i[1], i[2]], findfirst(x-> x == j, grid[i[1], i[2]])))
                elseif i[2] != index[2]
                    push!(grid[i[1], index[2]], splice!(grid[i[1], i[2]], findfirst(x-> x == j, grid[i[1], i[2]])))
                end
            end
        end
    end

    function force(p1, p2, σ, ϵ, r_cut)
        rmag = norm(p1.x - p2.x)
        r_vector = (p2.x - p1.x)/rmag ## Normalized vector
        if rmag <= r_cut
            s = (σ/rmag)^6
            return 24ϵ*(s)*(1-2s)*r_vector
        else
            return [0, 0]
        end
    end

    function compV_LC(grid, dt, nc)
        for i in Iterators.product(1:nc, 1:nc) |> collect
            for k in grid[i[1], i[2]]
                updateV(k, dt)
            end
        end
    end

    function updateV(p, dt)
        p.v = p.v + (dt*0.5/p.m)*(p.F + p.F_old)
    end

    function run(particles, dt, t_end, r_cut, σ, ϵ, L, nc)
        N, n = length(particles), Int(t_end/dt)
        grid = grid_init(particles, nc, L, r_cut)
        result_x, result_y, result_e = timeIntegration_LC(0, dt, t_end, grid, N, n, r_cut, σ, ϵ, L, nc)
        return result_x, result_y, result_e
    end

    function grid_init(particles, nc, L, r_cut)
        grid = Array{Any}(undef, nc, nc)
        for i in Iterators.product(1:nc, 1:nc) |> collect
            grid[i[1], i[2]] = Array{Any}(undef, 0)
        end
        for i in particles
            index = div.(i.x, r_cut, RoundUp) .|> Int
            push!(grid[index[1], index[2]], i)
        end
        return grid
    end
    
    ### Auxiliary functions ###

    function legal(x)
        return filter!(e->1<=e<=10, x-1:x+1 |> collect)
    end
    
    function neighbour_cells(x, y)
        return Iterators.product(legal(x)[1]:legal(x)[end], legal(y)[1]:legal(y)[end]) |> collect
    end

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
end
