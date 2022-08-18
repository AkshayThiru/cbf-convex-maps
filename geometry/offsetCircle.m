% For systems with three states: [x-position, y-position, heading]
classdef offsetCircle < geometry
    properties
        center;
        radius;
        offset;
    end
    
    methods
        function obj = offsetCircle(c, r, o)
            % c: Distance of the center of the circle [m].
            % r: Radius of the circle [m].
            % o: Offset angle of the center from the robot heading [rad].
            obj@geometry(2, 1);
            
            obj.center = c;
            obj.radius = r;
            obj.offset = o;
        end
        
        function vec = A(obj, x, z)
            vec = (z(1)-x(1)-obj.center*cos(x(3)+obj.offset))^2 + ...
                (z(2)-x(2)-obj.center*sin(x(3)+obj.offset))^2 - obj.radius^2;
        end
        
        function mat = dAdx(obj, x, z)
            mat = [2*(x(1)+obj.center*cos(x(3)+obj.offset)-z(1)) ...
                2*(x(2)+obj.center*sin(x(3)+obj.offset)-z(2)) ...
                2*obj.center*((z(1)-x(1))*sin(x(3)+obj.offset)-(z(2)-x(2))*cos(x(3)+obj.offset))];
        end
        
        function mat = dAdz(obj, x, z)
            mat = [2*(z(1)-x(1)-obj.center*cos(x(3)+obj.offset)) ...
                2*(z(2)-x(2)-obj.center*sin(x(3)+obj.offset))];
        end
        
        function mat = d2Adxz_L(obj, x, ~, L)
            mat = L(1)*[-2 0 2*obj.center*sin(x(3)+obj.offset); ...
                0 -2 -2*obj.center*cos(x(3)+obj.offset)];
        end
        
        function mat = d2Adz2_L(~, ~, ~, L)
            mat = L(1)*[2 0; 0 2];
        end
        
        function [points] = plot_outline(obj, x, ang_lim)
            N = 200;
            if nargin == 3
                N = round(N/2/pi* acos(cos(ang_lim(2)-ang_lim(1))));
                temp = linspace(ang_lim(1), ang_lim(2), N);
            elseif nargin == 2
                temp = linspace(0, 2*pi, N);
            end
            points = [obj.radius*cos(temp)+obj.center; obj.radius*sin(temp)];
            points = x(1:2)*ones(1,N) + ...
                [cos(x(3)+obj.offset) -sin(x(3)+obj.offset); sin(x(3)+obj.offset) cos(x(3)+obj.offset)]*points;
        end
    end
end
