classdef (Abstract) robot
    properties
        n % Dim. of x (state).
        m % Dim. of u (input).
        l % Dim. of z (position).
        r % Num of inequalities.
        
        x  % State.
        xf % Desired final state.
        
        A        % Handle for inequality: A(x, z): (r, 1).
        dAdx     % Handle: dAdx(x, z): (r, n).
        dAdz     % Handle: dAdz(x, z): (r, l).
        d2Adxz_L % Handle: d2Adxz_L(x, z, L): (l, n): \sum_{i=1}^r Li * d2Ai/dxz.
        d2Adz2_L % Handle: d2Adz2_L(x, z, L): (l, l): \sum_{i=1}^r Li * d2Ai/dz2.
        
        geo; % Geometry object.
    end
    
    methods
        function obj = robot(n, m)
            obj.n = n;
            obj.m = m;
            
            obj.x = zeros(n, 1);
        end
        
        function xdot = dyn(obj, x, u)
            xdot = obj.f(x) + obj.g(x)*u;
        end
        
        function obj = assign_geometry(obj, geo)
            obj.l = geo.l;
            obj.r = geo.r;
            
            obj.A        = @(x,z) geo.A(x,z);
            obj.dAdx     = @(x,z) geo.dAdx(x,z);
            obj.dAdz     = @(x,z) geo.dAdz(x,z);
            obj.d2Adxz_L = @(x,z,L) geo.d2Adxz_L(x,z,L);
            obj.d2Adz2_L = @(x,z,L) geo.d2Adz2_L(x,z,L);
            
            obj.geo = geo;
        end
    end
    
    methods (Abstract)
        f(obj, x) % Drift vector field: f(x): (n, 1).
        g(obj, x) % Input vector field: g(x): (n, m).
    end
end
