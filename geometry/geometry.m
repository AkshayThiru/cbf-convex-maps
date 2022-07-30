classdef geometry
    properties
        l % Dim. of z (position).
        r % Num of inequalities.
        
        A        % Handle for inequality: A(x, z): (r, 1).
        dAdx     % Handle: dAdx(x, z): (r, n).
        dAdz     % Handle: dAdz(x, z): (r, l).
        d2Adxz_L % Handle: d2Adxz_L(x, z, L): (l, n): \sum_{i=1}^r Li * d2Ai/dxz.
        d2Adz2_L % Handle: d2Adz2_L(x, z, L): (l, l): \sum_{i=1}^r Li * d2Ai/dz2.
    end
    
    methods
        function obj = geometry(l, r)
            obj.l = l;
            obj.r = r;
        end
    end
end
