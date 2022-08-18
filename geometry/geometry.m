classdef (Abstract) geometry
    properties
        l % Dim. of z (position).
        r % Num of inequalities.
    end
    
    methods
        function obj = geometry(l, r)
            obj.l = l;
            obj.r = r;
        end
    end
    
    methods (Abstract)
        A(obj, x, z) % Inequality constraints: A(x, z): (r, 1).
        dAdx(obj, x, z) % dAdx(x, z): (r, n).
        dAdz(obj, x, z) % dAdz(x, z): (r, l).
        d2Adxz_L(obj, x, z, L) % d2Adxz_L(x, z, L): (l, n): \sum_{i=1}^r Li * d2Ai/dxz.
        d2Adz2_L(obj, x, z, L) % d2Adz2_L(x, z, L): (l, l): \sum_{i=1}^r Li * d2Ai/dz2.
    end
end
