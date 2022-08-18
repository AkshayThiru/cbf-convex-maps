function [arr] = init_robots()
    %% Define robots.
    s = 3;
    
    switch s
        case 1
            arr = cell(1,2);

            % Robot 1.
            r = integrator();
            r = r.assign_geometry(nCircles(2, 0.5, 1, 0*pi/180));
            r.x  = [0; 1; 0];
            r.xf = [5; 0; 0];
            arr{1} = r;

            % Robot 2.
            r = unicycle();
            r = r.assign_geometry(nCircles(2, 0.5, 1, 90*pi/180));
            r.x  = [5; 0; pi];
            r.xf = [0; 0; 45.0*pi/180];
            arr{2} = r;
            
        case 2
            n = 5;
            R = 10;
            
            arr = cell(1, n);
            for i = 1:n
                r = integrator();
%                 r = r.assign_geometry(offsetCircle(0, 1, 0));
%                 r = r.assign_geometry(nCircles(3, 0.5, 1, 0));
                r = r.assign_geometry(highOrderEllipse(4, 1, 1, 0));
                r.x  = [R*cos(2*pi*i/n); R*sin(2*pi*i/n); 0];
                r.xf = [R*cos(2*pi*i/n+pi); R*sin(2*pi*i/n+pi); 0];
                arr{i} = r;
            end
        
        case 3
            n = 5;
            R = 9;
            
            arr = cell(1, n);
            % Degree 4 ellipse.
            arr{1} = integrator();
            arr{1} = arr{1}.assign_geometry(highOrderEllipse(4, 1.5, 1, 0));
            % 3 circle intersection.
            arr{2} = integrator();
            arr{2} = arr{2}.assign_geometry(nCircles(3, 1, 2, 0));
            % 2 circle intersection.
            arr{3} = integrator();
            arr{3} = arr{3}.assign_geometry(nCircles(2, 1, 2, 0));
            % Ellipse.
            arr{4} = unicycle();
            arr{4} = arr{4}.assign_geometry(highOrderEllipse(2, 2, 1, 0));
            % Degree 4 ellipse.
            arr{5} = unicycle();
            arr{5} = arr{1}.assign_geometry(highOrderEllipse(6, 2, 1, 0));
            % Offset Circle.
            if n >= 6
                arr{6} = unicycle();
                arr{6} = arr{5}.assign_geometry(offsetCircle(1.5, 1, 0));
            end
            
            for i = 1:n
                arr{i}.x  = [R*cos(2*pi*i/(n+0)); 1/2*R*sin(2*pi*i/(n+0)); 0];
                arr{i}.xi = arr{i}.x;
                arr{i}.xf = [R*cos(2*pi*i/(n+0)+pi); 1/2*R*sin(2*pi*i/(n+0)+pi); 0];
            end
    end
end

