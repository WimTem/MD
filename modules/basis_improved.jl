module basis_improved
    using LinearAlgebra

    function timeIntegration_basis(t::Real, dt::Real, t_end::Real, p, N::Int64, n::Int64, r_cut::Real, σ::Real, ϵ::Real, L::Real)
        result_x = zeros(n, N)
        result_y = zeros(n, N)
        result_e = zeros(n)
        result_u = zeros(n)
        count = 0
        for i = 1:n
            count += 1
            compX_basis(p, N, dt, L)
            compF_basis(p, N, r_cut, σ, ϵ)
            compV_basis(p, N, dt)
            for j = 1:N
                result_x[count, j] = p[j].x[1]
                result_y[count, j] = p[j].x[2]
                temp = compoutStatistic_basis(p, N, σ, ϵ, r_cut)
                result_e[count], result_u[count] = temp[1], temp[2]
            end
        end
        return result_x, result_y, result_e, result_u
    end

    function compoutStatistic_basis(p, N::Int64, σ, ϵ, r_cut)
        e, u = 0, 0
        for i = 1:N
            e += 0.5*sum(p[i].v)^2*p[i].m
        end
        for i = 1:N
            for j = 1:N
                if i < j
                    rmag = norm(p[i].x - p[j].x)
                    if rmag < r_cut^2
                        u += 4*ϵ*((σ/rmag)^12 - (σ/rmag)^6)
                    end
                end
            end
        end
        return e, u
    end

    function compX_basis(p, N::Int64, dt::Real, L::Real)
        for i = 1:N
            updateX(p[i], dt)
            ##Periodic boundaries
            for j = 1:2
                if p[i].x[j] <= 0
                    p[i].x[j] += L
                end
                if p[i].x[j] > L
                    p[i].x[j] -= L
                end
            end
        end
    end

    function updateX(p, dt)
        p.x += dt*p.v + (dt^2)*(1/2).*p.F
        p.F_old = p.F
    end

    function compF_basis(p, N::Int64, r_cut, σ, ϵ)
        for i = 1:N
            p[i].F = [0,0]
        end
        for i = 1:N
            for j = 1:N
                if i < j
                    f = force(p[i], p[j], r_cut, σ, ϵ)
                    p[i].F = f
                    p[j].F = -f
                end
            end
        end
    end

    function force(p1, p2, r_cut, σ, ϵ)
        rmag = norm(p1.x - p2.x)
        r_vector = (p1.x - p2.x)/rmag ## Normalized vector
        if rmag <= r_cut
            r2 = 1/rmag^2
            r6 = r2^3
            return 48*r2*r6*ϵ*(r6 - .5)*r_vector
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