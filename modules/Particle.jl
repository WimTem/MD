module Particle
    mutable struct Part
        m
        x
        v
        F
        F_old
        Part(m, x, v, F, F_old) = new(m, x, v, F, F_old)
    end
end