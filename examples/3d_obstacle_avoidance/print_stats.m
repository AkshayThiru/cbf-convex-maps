function [] = print_stats(logs)
    nrobots = logs.nrobots;
    nT = logs.nT;

    dist_solve_time = zeros(nrobots, nT - 2);
    qp_solve_time = logs.solve_time_seq{1}(2, :);
    for i = 1:nrobots
        dist_solve_time(i, :) = logs.solve_time_seq{i}(1, 2:end);
    end
    p_ = pctl(dist_solve_time(:), [50, 99]);
    disp(['Distance ODE solution time (s): ', ...
        'mean = ', num2str(mean(dist_solve_time(:))), ', ', ...
        'std = ', num2str(std(dist_solve_time(:))), ', ', ...
        'p50 = ', num2str(p_(1)), ', ', ...
        'p99 = ', num2str(p_(2)), ', ', ...
        ]);
    p_ = pctl(qp_solve_time(:), [50, 99]);
    disp(['QP solution time (s): ', ...
        'mean = ', num2str(mean(qp_solve_time(:))), ', ', ...
        'std = ', num2str(std(qp_solve_time(:))), ', ', ...
        'p50 = ', num2str(p_(1)), ', ', ...
        'p99 = ', num2str(p_(2)), ', ', ...
        ]);

    dist2_ode_seq = zeros(nrobots * (nrobots - 1) / 2, nT);
    dist2_opt_seq = zeros(nrobots * (nrobots - 1) / 2, nT);
    idx_ = 1;
    for i = 1:nrobots
        for j = i+1:nrobots
            dist2_ode_seq(idx_, :) = logs.dist2_seq{i, j};
            dist2_opt_seq(idx_, :) = logs.opt_dist2_seq{i, j};
            idx_ = idx_ + 1;
        end
    end
    dist_err_seq = abs(sqrt(dist2_opt_seq(:)) - sqrt(dist2_ode_seq(:)));
    dist_err = max(dist_err_seq);
    dist_rel_err = max(dist_err_seq ./ sqrt(dist2_opt_seq(:)));
    disp(['Distance Inf norm error (m): ' num2str(dist_err)]);
    disp(['Distance Inf relative error (): ' num2str(dist_rel_err)]);

    dist_opt_seq_ = sqrt(dist2_opt_seq(:, 1:end-1));
    dist_ode_seq_ = sqrt(dist2_ode_seq(:, 1:end-1));
    Ddist_err_seq = abs(logs.opt_Ddist2_seq(:) ./ (2*dist_opt_seq_(:)) - ...
        logs.Ddist2_seq(:) ./ (2*dist_ode_seq_(:)));
    p_ = pctl(Ddist_err_seq, [50, 99]);
    disp(['Dist derivative error (m/s): ', ...
        'mean = ', num2str(mean(Ddist_err_seq)), ', ', ...
        'std = ', num2str(std(Ddist_err_seq)), ', ', ...
        'p50 = ', num2str(p_(1)), ', ', ...
        'p99 = ', num2str(p_(2)), ', ', ...
        ]);
    Ddist_ode_seq_ = abs(logs.opt_Ddist2_seq(:) ./ (2*dist_opt_seq_(:)));
    disp(['Dist derivative Inf (m/s): ', num2str(max(Ddist_ode_seq_))]);

    zopt_ode_seq = zeros(nrobots * (nrobots - 1) / 2, 8, nT);
    zopt_opt_seq = zeros(nrobots * (nrobots - 1) / 2, 8, nT);
    yopt_ode_seq = zeros(nrobots * (nrobots - 1) / 2, 2, nT);
    yopt_opt_seq = zeros(nrobots * (nrobots - 1) / 2, 2, nT);
    idx_ = 1;
    for i = 1:nrobots
        for j = i+1:nrobots
            zopt_ode_seq(idx_, :, :) = logs.zopt_seq{i, j};
            zopt_opt_seq(idx_, :, :) = logs.opt_zopt_seq{i, j};
            yopt_ode_seq(idx_, :, :) = logs.yopt_seq{i, j};
            yopt_opt_seq(idx_, :, :) = logs.opt_yopt_seq{i, j};
            idx_ = idx_ + 1;
        end
    end
    zopt_err_seq = vecnorm(zopt_opt_seq - zopt_ode_seq, 2, 2);
    yopt_err_seq = vecnorm(yopt_opt_seq - yopt_ode_seq, 2, 2);
    zopt_opt_norm_seq = vecnorm(zopt_opt_seq, 2, 2);
    yopt_opt_norm_seq = vecnorm(yopt_opt_seq, 2, 2);
    zopt_rel_err_seq = zopt_err_seq(:) ./ zopt_opt_norm_seq(:);
    yopt_rel_err_seq = yopt_err_seq(:) ./ yopt_opt_norm_seq(:);
    zopt_rel_err = max(zopt_rel_err_seq);
    yopt_rel_err = max(yopt_rel_err_seq);
    disp(['Primal solution Inf relative error (): ' num2str(zopt_rel_err)]);
    disp(['Dual solution Inf relative error (): ' num2str(yopt_rel_err)]);

    kkt_err_seq_ = max(logs.kkt_err_seq, [], 1);
    p_ = pctl(kkt_err_seq_, [50, 99]);
    disp(['KKT error (): ', ...
        'mean = ', num2str(mean(kkt_err_seq_)), ', ', ...
        'std = ', num2str(std(kkt_err_seq_)), ', ', ...
        'p50 = ', num2str(p_(1)), ', ', ...
        'p99 = ', num2str(p_(2)), ', ', ...
        ]);
    disp(['KKT error Inf (): ', num2str(max(kkt_err_seq_))]);
end
