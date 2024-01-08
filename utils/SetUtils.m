classdef SetUtils
    properties (Constant)
        MIN_POLYTOPE_INTERIOR_RADIUS = 1e-4; % [m]
        MAXMIN_CONVEXSET_FUNC_VALUE  = -1e-3;
    end
    
    methods (Static)
        function is_solid = is_polytope_solid(A, b)
            P = Polyhedron(A, b);
            s = P.chebyCenter;
            is_solid = (s.r >= SetUtils.MIN_POLYTOPE_INTERIOR_RADIUS);
        end
        
        function is_bounded = is_polytope_bounded(A, b)
            P = Polyhedron(A, b);
            is_bounded = P.isBounded;
        end
        
        function is_solid = is_convexset_solid(nz, C, x_test)
            obj = @(z) min(C.A(x_test, z), 1);
            z0 = zeros(nz, 1);
            
            options = optimoptions('fminunc', 'Display', 'off', ...
                'Algorithm', 'quasi-newton');
            options.OptimalityTolerance = 1e-7;
            [~, fval] = fminunc(obj, z0, options);
            
            is_solid = (fval <= SetUtils.MAXMIN_CONVEXSET_FUNC_VALUE);
        end
    end
end
