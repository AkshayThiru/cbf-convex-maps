% For systems with three states: [x-position, y-position, heading]
classdef two_circles < geometry
    properties
        c1;
        c2;
        r1;
        r2;
        o;
    end
    
    methods
        function obj = two_circles(c1, c2, r1, r2, o)
            % ci: Offset distance of the center of circle i [m].
            % ri: Radius of circle i [m].
            % o : Offset angle of the geometry from the robot heading [rad].
            obj@geometry(2, 2);
            
            obj.c1 = c1;
            obj.c2 = c2;
            obj.r1 = r1;
            obj.r2 = r2;
            obj.o  = o;
            
            obj.A = @(x, z) [(z(1)-x(1)-c1*cos(x(3)+o))^2 + (z(2)-x(2)-c1*sin(x(3)+o))^2 - r1^2;
                (z(1)-x(1)+c2*cos(x(3)+o))^2 + (z(2)-x(2)+c2*sin(x(3)+o))^2 - r2^2];
            obj.dAdx = @(x, z) [2*(x(1)+c1*cos(x(3)+o)-z(1)) 2*(x(2)+c1*sin(x(3)+o)-z(2)) 2*c1*(z(1)-x(1))*sin(x(3)+o)-2*c1*(z(2)-x(2))*cos(x(3)+o);
                2*(x(1)-c2*cos(x(3)+o)-z(1)) 2*(x(2)-c2*sin(x(3)+o)-z(2)) -2*c2*(z(1)-x(1))*sin(x(3)+o)+2*c2*(z(2)-x(2))*cos(x(3)+o)];
            obj.dAdz = @(x, z) [2*(z(1)-x(1)-c1*cos(x(3)+o)) 2*(z(2)-x(2)-c1*sin(x(3)+o)); ...
                2*(z(1)-x(1)+c2*cos(x(3)+o)) 2*(z(2)-x(2)+c2*sin(x(3)+o))];
            obj.d2Adxz_L = @(x, z, L) L(1)*[-2 0 2*c1*sin(x(3)+o); 0 -2 -2*c1*cos(x(3)+o)] + ...
                L(2)*[-2 0 -2*c2*sin(x(3)+o); 0 -2 2*c2*cos(x(3)+o)];
            obj.d2Adz2_L = @(x, z, L) (L(1)+L(2))*[2 0; 0 2];
        end
        
        function [points] = plot_outline(obj, x)
            N = 100;
            x_int = 1/2 * ((obj.r2^2-obj.r1^2)/(obj.c1+obj.c2) + obj.c1-obj.c2);
            y_int = sqrt(abs(obj.r1^2 - (x_int-obj.c1)^2));
            ang_1 = atan2(y_int, x_int+obj.c2);
            temp = linspace(-ang_1, ang_1, N);
            points = [obj.r2*cos(temp)-obj.c2; obj.r2*sin(temp)];
            ang_2 = atan2(y_int, x_int-obj.c1);
            temp = linspace(ang_2, 2*pi-ang_2, N);
            points = [points [obj.r1*cos(temp)+obj.c1; obj.r1*sin(temp)]];
            points = x(1:2)*ones(1,2*N) + ...
                [cos(x(3)+obj.o) -sin(x(3)+obj.o); sin(x(3)+obj.o) cos(x(3)+obj.o)]*points;
        end
    end
end
