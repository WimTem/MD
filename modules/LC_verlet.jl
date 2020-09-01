 module LC_verlet
    using LinearAlgebra, Random, Test

    function timeIntegration_LC(t::Real, dt::Real, t_end::Real, grid, N::Int64, n::Int64, r_cut, σ, ϵ, L, nc)
        result_x = zeros(n, N)
        result_y = zeros(n, N)
        result_e = zeros(n)
        result_u = zeros(n)
        count = 0
        for i = 1:n
            count += 1
            compX_LC(grid, dt, r_cut, L, nc)
            u = compF_LC(grid, σ, ϵ, r_cut, nc)
            compV_LC(grid, dt, nc)
            z = 0
            result_u[count] += u
            for i in Iterators.product(1:nc, 1:nc) |> collect
                if length(grid[i[1], i[2]]) > 0 && count <= n
                    for k in grid[i[1], i[2]]
                        z += 1
                        result_x[count, z] = k.x[1]
                        result_y[count, z] = k.x[2]
                        result_e[count] += compoutKineticEnergy_LC(k)
                    end
                end
            end
        end
        return result_x, result_y, result_e, result_u
    end

    function compoutKineticEnergy_LC(k)
        return .5*(k.v[1]^2 + k.v[2]^2)/k.m
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
        p.x = p.x + dt*p.v + dt*dt*(p.F)*(1/p.m)
        p.F_old = p.F
    end

    function compF_LC(grid, σ, ϵ, r_cut, nc)
        u = 0
        for i in Iterators.product(1:nc, 1:nc) |> collect
            for p1 in grid[i[1], i[2]]
                for j in neighbour_cells(i[1], i[2], nc)
                    for p2 in grid[j[1], j[2]]
                        if p1 != p2
                            p1.F = force(p1, p2, σ, ϵ, r_cut)
                            r = norm(p1.x-p2.x)
                            u += 4*ϵ*((σ/r)^12 - (σ/r)^6) - 4*ϵ*((σ/r_cut)^12 - (σ/r_cut)^6)
                        end
                    end
                end
            end
        end
        return u
    end

    function moveParticles_LC(grid, nc, L, r_cut)
        for i in Iterators.product(1:nc, 1:nc) |> collect
            for j in grid[i[1], i[2]]
            # Periodic boundaries
                for k = 1:2
                    if j.x[k] <= 0
                        j.x[k] = (j.x[k] % L) + L
                    elseif j.x[k] > L
                        j.x[k] = j.x[k] % L 
                    end
                end
            end
            ## Move back to correct cell
            for j in grid[i[1], i[2]]
                index = (div.(j.x, r_cut, RoundUp)) .|> Int
                for k = 1:2
                    if index[k] == 0
                        index[k] = 1
                    elseif index[k] == (nc+1)
                        index[k] = nc
                    end
                end
                if i[1] != index[1]
                    push!(grid[index[1], i[2]], splice!(grid[i[1], i[2]], findfirst(x-> x == j, grid[i[1], i[2]])))
                elseif i[2] != index[2]
                    push!(grid[i[1], index[2]], splice!(grid[i[1], i[2]], findfirst(x-> x == j, grid[i[1], i[2]])))
                end
            end
        end
    end

    function force(p1, p2, r_cut, σ, ϵ)
        rmag = norm(p1.x - p2.x)
        r_vector = (rmag^-1).*(p1.x - p2.x) ## Normalized vector
        if rmag <= r_cut
            r6 = (σ/rmag)^6
            return 24*ϵ/rmag*(2*r6^2 - r6).*r_vector
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
        p.v = p.v + 0.5* (p.F)*dt/p.m
    end

    function run(particles, dt, t_end, r_cut, σ, ϵ, L, nc)
        N, n = length(particles), Int(t_end/dt)
        grid = grid_init(particles, nc, L, r_cut)
        result_x, result_y, result_e, result_u = timeIntegration_LC(0, dt, t_end, grid, N, n, r_cut, σ, ϵ, L, nc)
        return result_x, result_y, result_e, result_u
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
    function legal(x, nc)
        return filter!(e->1<=e<=nc, x-1:x+1 |> collect)
    end
    
    function neighbour_cells(x, y, nc)
        return Iterators.product(legal(x, nc)[1]:legal(x, nc)[end], legal(y, nc)[1]:legal(y, nc)[end]) |> collect
    end
        
end