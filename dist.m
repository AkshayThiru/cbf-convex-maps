function [hij, z, L, J] = dist(ri, rj, z_init)
    eps = 1e-6; % Margin for index calculations.
    
    %% Minimum distance between strictly convex sets.
    l = ri.l; % = rj.l
    
    % Objective function.
    obj  = @(z) objective(z, l);
    cons = @(z) nonlcon(z, ri, rj);
    
    % Constraint bounds.
    options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    options.SpecifyConstraintGradient = true;
    options.SpecifyObjectiveGradient  = true;
    options.ConstraintTolerance       = 1e-7;
    options.OptimalityTolerance       = 1e-7;
%     options.HessianFcn                = @(z, lambda) hessianfunc(z, lambda, ri, rj);
%     options.CheckGradients            = true;
    
    % Initialize the primal variables.
    if exist('z_init', 'var')
        z0 = z_init;
    else
        z0 = zeros(2*l, 1);
    end
    
    % Solve the minimum distance problem.
    % tic;
    [z_opt, hij, ~,~, lambda, ~,~] = fmincon(obj, z0, [],[],[],[],[],[], cons, options);
    % toc
    
    % Extract outputs.
    z.i = z_opt(1:l);
    z.j = z_opt(l+1:2*l);
    L.i = lambda.ineqnonlin(1:ri.r);
    L.j = lambda.ineqnonlin(ri.r+1:ri.r+rj.r);
    
    J.J_0c = find([ri.A(ri.x, z.i); rj.A(rj.x, z.j)] < -eps);
    J.J_1  = setdiff(find([L.i; L.j] > eps), J.J_0c);
    J.J_2e = setdiff((1:ri.r+rj.r)', union(J.J_0c, J.J_1));
end

function [fun, grad] = objective(z, l)
    fun = (z(1:l)-z(l+1:2*l))' * (z(1:l)-z(l+1:2*l));
    grad = 2*[z(1:l)-z(l+1:2*l); z(l+1:2*l)-z(1:l)];
end

function [c, ceq, GC, GCeq] = nonlcon(z, ri, rj)
    l = ri.l;
    
    c = [ri.A(ri.x, z(1:l)); rj.A(rj.x, z(l+1:2*l))]';
    ceq = [];
    
    GC = [ri.dAdz(ri.x, z(1:l)) zeros(ri.r, l); zeros(rj.r, l) rj.dAdz(rj.x, z(l+1:2*l))]';
    GCeq = [];
end

function H = hessianfunc(z, lambda, ri, rj)
    l = ri.r;
    v = lambda.ineqnonlin;
    
    H = [2*eye(l)+ri.d2Adz2_L(ri.x, z(1:l), v(1:ri.r)) -2*eye(l);
        -2*eye(l) 2*eye(l)+rj.d2Adz2_L(rj.x, z(l+1:2*l), v(ri.r+1:ri.r+rj.r))];
end
