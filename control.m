function [u, prev, log, debug] = control(arr, prev, const)
    global safety
    
    A_cbf = zeros(const.N, const.M);
    b_cbf = zeros(const.N, 1);
    
    log = struct('h', [], 't_opt', []);
    debug = [];
    
    for i = 1:length(arr)
        for j = i+1:length(arr)
            k = ((2*length(arr)-i)*(i-1))/2 + j-i;
            
            [Aij, bij, prev{k}, log_ij, ~] = input_constraint_ij(arr{i}, arr{j}, prev{k});
            A_cbf(k, const.m(i):const.m(i+1)-1) = Aij(1:arr{i}.m);
            A_cbf(k, const.m(j):const.m(j+1)-1) = Aij(1+arr{i}.m:end);
            b_cbf(k) = bij;
            
            % Logging.
            logn = fieldnames(log);
            for w = 1:numel(logn)
                log.(logn{w})(:,k) = log_ij.(logn{w});
            end
        end
    end
    
    % Nominal controller.
    u_nom = zeros(const.M, 1);
    for i = 1:length(arr)
        u_nom(const.m(i):const.m(i+1)-1) = arr{i}.nominal_control(arr{i}.x);
    end
    
    H = 2 * eye(const.M);
    f = -u_nom;
    
    options = optimoptions('quadprog','Display','off');
    if ~safety
        % Nominal controller.
        tic;
        sol = quadprog(H, f, [],[],[],[],[],[],[], options);
        log.t_opt(:,const.N+1) = toc;
    else
        % Safety-critical controller.
        tic
        sol = quadprog(H, f, A_cbf, b_cbf, [],[],[],[],[], options);
        log.t_opt(:,const.N+1) = toc;
    end
    
    u = cell(length(arr), 1);
    for i = 1:length(arr)
        u{i} = sol(const.m(i):const.m(i+1)-1);
    end
end

function [A_cbf, b_cbf, prev_zij, log, debug] = input_constraint_ij(ri, rj, prev_zij)
    a_cbf = 1; % ECBF rate.
    e     = 0.1; % Safety margin [m].
    M_L   = 1e3; % Bound on dL. |dL|.|dA| ~ |dh/dt|.
    
    l = ri.l;
    
    % Compute minimum distances.
    [hij, z, L, J, t_opt] = dist(ri, rj, prev_zij);
    zi = z.i; zj = z.j;
    Li = L.i; Lj = L.j;
    prev_zij = [zi; zj];
    
    % Compute matrices.
    hess = [2*eye(l)+ri.d2Adz2_L(ri.x, zi, Li) -2*eye(l);
        -2*eye(l) 2*eye(l)+rj.d2Adz2_L(rj.x, zj, Lj)];
    hess_inv = hess^(-1);
    dA = [ri.dAdz(ri.x, zi)' zeros(l, rj.r);
        zeros(l, ri.r) rj.dAdz(rj.x, zj)'];
    mat_u_eq = [ri.d2Adxz_L(ri.x, zi, Li)*ri.g(ri.x) zeros(l, rj.m);
        zeros(l, ri.m) rj.d2Adxz_L(rj.x, zj, Lj)*rj.g(rj.x)];
    vec_c_eq = [ri.d2Adxz_L(ri.x, zi, Li)*ri.f(ri.x);
        rj.d2Adxz_L(rj.x, zj, Lj)*rj.f(rj.x)];
    vec_z_in = -[2*(zi-zj)'+Li'*ri.dAdz(ri.x, zi) 2*(zj-zi)'+Lj'*rj.dAdz(rj.x, zj)];
    vec_L_in = -[ri.A(ri.x, zi)' rj.A(rj.x, zj)'];
    vec_u_in = -[Li'*ri.dAdx(ri.x, zi)*ri.g(ri.x) Lj'*rj.dAdx(rj.x, zj)*rj.g(rj.x)];
    c_in = a_cbf*(hij-e) + Li'*ri.dAdx(ri.x, zi)*ri.f(ri.x) + Lj'*rj.dAdx(rj.x, zj)*rj.f(rj.x);
    
    % Fourier-Motzkin elimination.
    vec_L_cbf = vec_L_in - vec_z_in*hess_inv*dA;
    vec_u_cbf = vec_u_in - vec_z_in*hess_inv*mat_u_eq;
    c_fme = sum(abs(vec_L_cbf(J.J_1))) * M_L + ...
        sum(-vec_L_cbf( intersect(find(vec_L_cbf < 0), J.J_2e) ))*M_L;
    
    % Compute input constraint.
    A_cbf = vec_u_cbf;
    b_cbf = c_in + vec_z_in*hess_inv*vec_c_eq + c_fme;
    
    % Logging
    log = struct('h', hij, 't_opt', t_opt);
    debug = [];
end
