 module LC_verlet
    using LinearAlgebra, Random

    function timeIntegration_LC(t::Real, dt::Real, t_end::Real, grid, N::Int64, n::Int64, r_cut, σ, ϵ, L, nc)
        result_x = zeros(n, N)
        result_y = zeros(n, N)
        count = 0
        while (t < t_end)
            count += 1
            t += dt
            z = 0
            compX_LC(grid, dt, r_cut, L, nc)
            compF_LC(grid, σ, ϵ, r_cut, nc)
            compV_LC(grid, dt, nc)
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

    function compX_LC(grid, dt, r_cut, L, nc)
        for i = 1:nc
            for j = 1:nc
                if length(grid[i,j]) > 0
                    for k in grid[i, j]
                        updateX(k, dt)
                    end
                end
            end
        end
        moveParticles_LC(grid, nc, L, r_cut)
    end

    function updateX(p, dt)
        a = dt*0.5/p.m
        p.x = p.x + dt*(p.v + a*p.F)
        p.F_old = p.F
    end

    function compV_LC(grid, dt, nc)
        for i = 1:nc
            for j = 1:nc
                for k in grid[i, j]
                    updateV(k, dt)
                end
            end
        end
    end

    function updateV(p, dt)
        a = dt*0.5/p.m
        p.v = p.v + a*(p.F + p.F_old)
    end

    function compF_LC(grid, σ, ϵ, r_cut, nc)
        for i in Iterators.product(1:nc, 1:nc) |> collect
            for p1 in grid[i[1], i[2]]
                for j in neighbour_cells(i[1], i[2])
                    for p2 in grid[j[1], j[2]]
                        force(p1, p2, σ, ϵ, r_cut)
                    end
                end
            end
        end
        #Loop over cells
        # for i = 1:nc
        #     for j = 1:nc
        #         #Loop over elements in cell
        #         if length(grid[i, j]) > 0
        #             #Loop over neighbourhood cells of each element
        #             for el1 in grid[i, j]
        #                 for k in i-1:i+1
        #                     for l in j-1:j+1
        #                         if k > 0 && k <= nc
        #                             if l > 0 && l <= nc
        #                                 for el2 in grid[k, l]
        #                                     if el1 != el2
        #                                         force(el1, el2, σ, ϵ, r_cut)
        #                                     end
        #                                 end
        #                             end
        #                         end
        #                     end
        #                 end
        #             end
        #         end
        #     end
        # end
    end

    function moveParticles_LC(grid, nc, L, r_cut)
        for i = 1:nc
            for j = 1:nc
                ### Reflecting boundaries
                if length(grid[i, j]) > 0
                    for k in grid[i, j]
                        for m = 1:2
                            if k.x[m] <= 0
                                k.v[m] = -k.v[m]
                                k.x[m] = 1e-3 #Small positive required, see line 147
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
    end

    function force(p1, p2, σ, ϵ, r_cut)
        rmag = norm(p1.x - p2.x)
        r_vector = (p2.x - p1.x)/rmag ## Normalized vector
        if rmag <= r_cut
            s = (σ/rmag)^6
            p1.F = 24ϵ*(s)*(1-2s)*r_vector
        end
    end

    function run(particles, dt, t_end, r_cut, σ, ϵ, L, nc)
        N = length(particles)
        grid = grid_init(particles, nc, L, r_cut)
        n = Int(t_end/dt)
        result_x, result_y = timeIntegration_LC(0, dt, t_end, grid, N, n, r_cut, σ, ϵ, L, nc)
        return result_x, result_y
    end

    function grid_init(particles, nc, L, r_cut)
        grid = Array{Any}(undef, nc, nc)
        for i = 1:nc
            for j = 1:nc
                grid[i, j] = Array{Any}(undef, 0)
            end
        end
        for i in particles
            push!(grid[1, 1], i)
        end
        moveParticles_LC(grid, nc, L, r_cut)
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
