classdef integrator < robot
    properties
        kv = 0.5;
    end
    
    methods
        function obj = integrator()
            obj@robot(3, 3);
        end
        
        function fx = f(~, ~)
            fx = [0; 0; 0];
        end
        
        function gx = g(~, ~)
            gx = eye(3);
        end
        
        function u_nom = nominal_control(obj, x)
            u_nom = [obj.kv*(obj.xf(1)-x(1)); obj.kv*(obj.xf(2)-x(2)); 0];
        end
    end
end
