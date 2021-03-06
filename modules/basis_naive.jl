module basis_naive
    using LinearAlgebra

    function timeIntegration_basis(t::Real, dt::Real, t_end::Real, p, N::Int64, n::Int64, r_cut::Real, σ::Real, ϵ::Real, L::Real)
        result_x = zeros(n, N)
        result_y = zeros(n, N)
        result_e = zeros(n)
        count = 0
        for i = 1:n
            count += 1
            compX_basis(p, N, dt, L)
            compF_basis(p, N, r_cut, σ, ϵ)
            compV_basis(p, N, dt)
            for j = 1:N
                result_x[count, j] = p[j].x[1]
                result_y[count, j] = p[j].x[2]
                result_e[count] = compoutStatistic_basis(p, N)
            end
        end
        return result_x, result_y, result_e
    end

    function compoutStatistic_basis(p, N::Int64)
        e = 0
        for i = 1:N
            v = 0
            for j = 1:2
                v += (p[i].v[j])^2
            end
            e += 0.5*p[i].m*v
        end
        return e
    end

    function compX_basis(p, N::Int64, dt::Real, L::Real)
        for i = 1:N
            updateX(p[i], dt)
            ##Reflecting boundaries
            for j = 1:2
                if p[i].x[j] <= 0
                    p[i].v[j] = -p[i].v[j]
                    p[i].x[j] = 1e-6
                end
                if p[i].x[j] > L
                    p[i].v[j] = -p[i].v[j]
                    p[i].x[j] = L
                end
            end
        end
    end

    function updateX(p, dt)
        p.x = p.x + dt*(p.v + (dt*0.5/p.m)*p.F)
        p.F_old = p.F
    end

    function compF_basis(p, N::Int64, r_cut, σ, ϵ)
        for i = 1:N
            p[i].F = [0,0]
        end
        for i = 1:N
            for j = 1:N
                if i != j
                    p[i].F = p[i]. F + force(p[i], p[j], r_cut, σ, ϵ)
                end
            end
        end
    end

    function force(p1, p2, r_cut, σ, ϵ)
        rmag = norm(p1.x - p2.x)
        r_vector = (p2.x - p1.x)/rmag ## Normalized vector
        if rmag <= r_cut
            s = (σ/rmag)^6
            return 24ϵ*(s)*(1-2s)*r_vector
        else
            return [0, 0]
        end
    end

    function compV_basis(p, N::Int64, dt::Real)
        for i = 1:N
            updateV(p[i], dt)
        end
    end

    function updateV(p, dt)
        p.v = p.v + (dt*0.5/p.m)*(p.F + p.F_old)
    end

    function run(particles, dt, t_end, r_cut, σ, ϵ, L)
        n = Int(t_end/dt)
        N = length(particles)
        return timeIntegration_basis(0, dt, t_end, particles, N, n, r_cut, σ, ϵ, L)
    end
end