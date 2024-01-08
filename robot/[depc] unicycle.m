classdef unicycle < robot
    properties
        kr = 0.5;
        ka = 2;
        h_theta = 1;
        set_final_head = false;
    end
    
    methods
        function obj = unicycle()
            obj@robot(3, 2);
        end
        
        function fx = f(~, ~)
            fx = [0; 0; 0];
        end
        
        function gx = g(~, x)
            gx = [cos(x(3)) 0; sin(x(3)) 0; 0 1];
        end
        
        function u_nom = nominal_control(obj, x)
            rad = sqrt((obj.xf(1:2)-x(1:2))'*(obj.xf(1:2)-x(1:2)));
            [alpha, theta] = unicycle.ref_angle([cos(x(3)); sin(x(3))], obj.xf(1:2)-x(1:2));
            if obj.set_final_head
                if abs(alpha) < 1e-6
                    w_nom = obj.kr*obj.h_theta*sin(obj.xf(3)-theta);
                else
                    w_nom = -obj.ka*alpha + obj.kr*sin(alpha)*cos(alpha)* ...
                        (alpha+obj.h_theta*sin(obj.xf(3)-theta))/alpha;
                end
            else
                w_nom = -obj.ka*alpha + obj.kr*sin(alpha)*cos(alpha);
            end
            u_nom = [obj.kr*rad*cos(alpha); w_nom];
        end
    end
    
    methods (Static)
        function [ang, theta] = ref_angle(v, v_ref)
            theta = -atan2(v_ref(2), v_ref(1));
            Rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            v_trans = Rot * v;
            ang = atan2(v_trans(2), v_trans(1));
        end
    end
end
