function [dist2_new, z_opt_new, y_opt_new, out] = ...
    minimum_distance_step(dt, x1, dx1, C1, x2, dx2, C2, z_opt, y_opt, alg)
    % Performs a single integration step.
    % Inputs:
    %   dt: timestep for integration.
    %   x1, x2: Current parameters for C1 and C2.
    %   dx1, dx2: Rate of change of parameters x1 and x2.
    %   C1, C2: Convex sets, at least one of them strictly convex.
    %   z_opt: Optimal primal solution for x=(x1, x2).
    %   y_opt: Optimal dual solution for x=(x1, x2).
    %   alg: either 'euler' (default) or 'ode15s'.
    % Outputs:
    %   dist2_new: dist2 after dt time.
    %   z_opt_new: z_opt after dt time.
    %   y_opt_new: y_opt after dt time.
    nz = C1.nz; % = C2.nz;

    if nargin < 10
        alg = 'euler';
    end
    
    if strcmp(alg, 'euler')
        [dz_opt, dy_opt, out] = ...
            minimum_distance_ode(x1, dx1, C1, x2, dx2, C2, z_opt, y_opt);
        z_opt_new = z_opt + dt * dz_opt;
        y_opt_new = y_opt + dt * dy_opt;
    else
        error('''ode15s'' method is not implemented');
    end

    dist2_new = sum((z_opt(1:nz) - z_opt(1+nz:end)).^2);
end
