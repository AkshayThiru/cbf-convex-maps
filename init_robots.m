function [arr] = init_robots()
    %% Define robots.
    s = 2;
    
    switch s
        case 1
            arr = cell(1,2);

            % Robot 1.
            r = integrator();
            r = r.assign_geometry(two_circles(0.5, 0.5, 1.0, 1.0, 0.0*pi/180));
            r.x  = [0; 1; 0];
            r.xf = [5.0; 0.0; 0.0];
            arr{1} = r;

            % Robot 2.
            r = unicycle();
            r = r.assign_geometry(two_circles(0.5, 0.5, 1.0, 1.0, 90.0*pi/180));
            r.x  = [5; 0; pi];
            r.xf = [0.0; 0.0; 45.0*pi/180];
            arr{2} = r;
            
        case 2
            n = 10;
            R = 10;
            
            arr = cell(1, n);
            for i = 1:n
                r = integrator();
                r = r.assign_geometry(two_circles(0.5, 0.5, 1.0, 1.0, 90.0*i*pi/180));
%                 r = r.assign_geometry(offset_circle(0.0, 1.0, 0.0));
                r.x  = [R*cos(2*pi*i/n); R*sin(2*pi*i/n); 0.0];
                r.xf = [R*cos(2*pi*i/n+pi); R*sin(2*pi*i/n+pi); 0.0];
                arr{i} = r;
            end
    end
end

