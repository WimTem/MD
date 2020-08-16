module Particle
    mutable struct p
        m
        x
        v
        F
        F_old
        tag
        function p(m, x, v, F, F_old) 
            return new(m, x, v, F, F_old, 0)
        end
    end
end