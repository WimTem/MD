module Particle
    mutable struct p
        m
        x
        v
        F
        F_old
        p(m, x, v, F, F_old) = new(m, x, v, F, F_old)
    end
end