% For systems with three states: [x-position, y-position, heading]
classdef highOrderEllipse < geometry
    properties
        p;
        major;
        minor;
        offset;
        
        funcs;
    end
    
    methods
        function obj = highOrderEllipse(p, a, b, o)
            % p: Degree of curve (must be even).
            % a: Semi-major axis length [m].
            % b: Semi-minor axis length [m].
            % o: Offset angle of the center from the robot heading [rad].
            obj@geometry(2, 1);
            
            obj.p = p;
            obj.major  = a;
            obj.minor  = b;
            obj.offset = o;
            
            obj.funcs = obj.generate_functions(p, a, b, o);
        end
        
        function funcs = generate_functions(~, p, a, b, o)
            funcs = cell(5, 1);
            
            syms x1 x2 x3 real
            syms z1 z2 real
            
            tht = x3 + o;
            z_trans = [cos(tht) sin(tht); -sin(tht) cos(tht)]*[z1-x1; z2-x2];
            A = (z_trans(1)/a)^p + (z_trans(2)/b)^p - 1;
            funcs{1} = matlabFunction(simplify(A), 'Vars', [x1,x2,x3,z1,z2]);
            
            dAdx = [diff(A, x1) diff(A, x2) diff(A, x3)];
            funcs{2} = matlabFunction(simplify(dAdx), 'Vars', [x1,x2,x3,z1,z2]);
            
            dAdz = [diff(A, z1) diff(A, z2)];
            funcs{3} = matlabFunction(simplify(dAdz), 'Vars', [x1,x2,x3,z1,z2]);
            
            d2Adxz = [diff(dAdx, z1); diff(dAdx, z2)];
            funcs{4} = matlabFunction(simplify(d2Adxz), 'Vars', [x1,x2,x3,z1,z2]);
            
            d2Adz2 = [diff(dAdz, z1); diff(dAdz, z2)];
            funcs{5} = matlabFunction(simplify(d2Adz2), 'Vars', [x1,x2,x3,z1,z2]);
        end
        
        function vec = A(obj, x, z)
            vec = obj.funcs{1}(x(1), x(2), x(3), z(1), z(2));
        end
        
        function mat = dAdx(obj, x, z)
            mat = obj.funcs{2}(x(1), x(2), x(3), z(1), z(2));
        end
        
        function mat = dAdz(obj, x, z)
            mat = obj.funcs{3}(x(1), x(2), x(3), z(1), z(2));
        end
        
        function mat = d2Adxz_L(obj, x, z, L)
            mat = obj.funcs{4}(x(1), x(2), x(3), z(1), z(2)) * L(1);
        end
        
        function mat = d2Adz2_L(obj, x, z, L)
            mat = obj.funcs{5}(x(1), x(2), x(3), z(1), z(2)) * L(1);
        end
        
        function [points] = plot_outline(obj, x)
            N = 200;
            temp = linspace(0, 2*pi, N);
            tht = x(3)+obj.offset;
            Rot = [cos(tht) -sin(tht); sin(tht) cos(tht)];
            stht = sin(temp);
            ctht = cos(temp);
            points = [obj.major*(abs(stht)).^(2/obj.p).*sign(stht); ...
                obj.minor*(abs(ctht)).^(2/obj.p).*sign(ctht)];
            points = x(1:2)*ones(1,N) + Rot*points;
        end
    end
end
