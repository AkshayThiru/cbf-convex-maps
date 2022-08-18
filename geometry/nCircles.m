% For systems with three states: [x-position, y-position, heading]
classdef nCircles < geometry
    properties
        num;
        center;
        radius;
        offset;
        
        circles;
    end
    
    methods
        function obj = nCircles(num, c, r, o)
            % num: Number of intersecting circles.
            % c  : Offset distance of the center of circles [m].
            % r  : Radius of circle [m].
            % o  : Offset angle of the geometry from the robot heading [rad].
            obj@geometry(2, num);
            
            obj.num    = num;
            obj.center = c;
            obj.radius = r;
            obj.offset = o;
            
            obj.circles = cell(num, 1);
            for i = 1:num
                ang = 2*pi/num*(i-1) + o;
                obj.circles{i} = offsetCircle(c, r, ang);
            end
        end
        
        function vec = A(obj, x, z)
            vec = zeros(obj.num, 1);
            for i = 1:obj.num
                vec(i) = obj.circles{i}.A(x, z);
            end
        end
        
        function mat = dAdx(obj, x, z)
            mat = zeros(obj.num, 3);
            for i = 1:obj.num
                mat(i,:) = obj.circles{i}.dAdx(x, z);
            end
        end
        
        function mat = dAdz(obj, x, z)
            mat = zeros(obj.num, 2);
            for i = 1:obj.num
                mat(i,:) = obj.circles{i}.dAdz(x, z);
            end
        end
        
        function mat = d2Adxz_L(obj, x, z, L)
            mat = zeros(2, 3);
            for i = 1:obj.num
                mat = mat + obj.circles{i}.d2Adxz_L(x, z, L(i));
            end
        end
        
        function mat = d2Adz2_L(obj, x, z, L)
            mat = zeros(2, 2);
            for i = 1:obj.num
                mat = mat + obj.circles{i}.d2Adz2_L(x, z, L(i));
            end
        end
        
        function [points] = plot_outline(obj, x)
            points = [];
            for i = 1:obj.num
                tht = pi/obj.num;
                chord_len = sqrt(obj.radius^2-obj.center^2*sin(tht)^2)-obj.center*cos(tht);
                if obj.num == 1
                    ang = pi;
                else
                    ang = asin(chord_len*sin(tht)/obj.radius);
                end
                ang_lim = obj.offset+pi + [-ang ang];
                points = [points obj.circles{i}.plot_outline(x, ang_lim)]; %#ok<AGROW>
            end
        end
    end
end
