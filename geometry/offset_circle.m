% For systems with three states: [x-position, y-position, heading]
classdef offset_circle < geometry
    properties
        c1;
        r1;
        o;
    end
    
    methods
        function obj = offset_circle(c1, r1, o)
            % c1: Distance of the center of the circle [m].
            % r1: Radius of the circle [m].
            % o : Offset angle of the center from the robot heading [rad].
            obj@geometry(2, 1);
            
            obj.c1 = c1;
            obj.r1 = r1;
            obj.o  = o;
            
            obj.A = @(x, z) [(z(1)-x(1)-c1*cos(x(3)+o))^2 + (z(2)-x(2)-c1*sin(x(3)+o))^2 - r1^2];
            obj.dAdx = @(x, z) [2*(x(1)+c1*cos(x(3)+o)-z(1)) 2*(x(2)+c1*sin(x(3)+o)-z(2)) 2*c1*(z(1)-x(1))*sin(x(3)+o)-2*c1*(z(2)-x(2))*cos(x(3)+o)];
            obj.dAdz = @(x, z) [2*(z(1)-x(1)-c1*cos(x(3)+o)) 2*(z(2)-x(2)-c1*sin(x(3)+o))];
            obj.d2Adxz_L = @(x, z, L) L(1)*[-2 0 2*c1*sin(x(3)+o); 0 -2 -2*c1*cos(x(3)+o)];
            obj.d2Adz2_L = @(x, z, L) L(1)*[2 0; 0 2];
        end
        
        function [points] = plot_outline(obj, x)
            N = 200;
            temp = linspace(0, 2*pi, N);
            points = [obj.r1*cos(temp)+obj.c1; obj.r1*sin(temp)];
            points = x(1:2)*ones(1,N) + ...
                [cos(x(3)+obj.o) -sin(x(3)+obj.o); sin(x(3)+obj.o) cos(x(3)+obj.o)]*points;
        end
    end
end
