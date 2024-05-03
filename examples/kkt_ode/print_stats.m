function [] = print_stats(out)
    disp(['Avg. opt solution time (s): ' num2str(mean(out.t_opt_seq))]);
    p_ = pctl(out.t_ode_seq(2:end), [50, 99]);
    disp(['ODE solution time (s): ', ...
        'mean = ', num2str(mean(out.t_ode_seq)), ', ', ...
        'std = ', num2str(std(out.t_ode_seq)), ', ', ...
        'p50 = ', num2str(p_(1)), ', ', ...
        'p99 = ', num2str(p_(2)), ', ', ...
        % 'inf = ', num2str(max(out.t_ode_seq(2:end))),
        ]);

    dist_err_seq = abs(sqrt(out.dist2_opt_seq) - sqrt(out.dist2_ode_seq));
    dist_err = max(dist_err_seq);
    dist_rel_err = max(dist_err_seq ./ sqrt(out.dist2_opt_seq));
    disp(['Distance Inf norm error (m): ' num2str(dist_err)]);
    disp(['Distance Inf relative error (): ' num2str(dist_rel_err)]);

    Ddist_err = abs(out.Ddist2_opt_seq ./ (2*sqrt(out.dist2_opt_seq(1:end-1))) - ...
        out.Ddist2_analytic_seq ./ (2*sqrt(out.dist2_ode_seq(1:end-1))));
    p_ = pctl(Ddist_err, [50, 99]);
    disp(['Dist derivative error (m/s): ', ...
        'mean = ', num2str(mean(Ddist_err)), ', ', ...
        'std = ', num2str(std(Ddist_err)), ', ', ...
        'p50 = ', num2str(p_(1)), ', ', ...
        'p99 = ', num2str(p_(2)), ', ', ...
        % 'inf = ', num2str(max(Ddist_err)),
        ]);
    Ddist_ode = abs(out.Ddist2_opt_seq ./ ...
        (2*sqrt(out.dist2_opt_seq(1:end-1))));
    disp(['Dist derivative Inf (m/s): ', num2str(max(Ddist_ode))]);

    zopt_err_seq = vecnorm(out.zopt_opt_seq - out.zopt_ode_seq, 2);
    zopt_rel_err_seq = zopt_err_seq ./ vecnorm(out.zopt_opt_seq, 2);
    zopt_rel_err = max(zopt_rel_err_seq);
    disp(['Primal solution Inf relative error (): ' num2str(zopt_rel_err)]);

    yopt_err_seq = vecnorm(out.yopt_opt_seq - out.yopt_ode_seq, 2);
    yopt_rel_err_seq = yopt_err_seq ./ vecnorm(out.yopt_opt_seq, 2);
    yopt_rel_err = max(yopt_rel_err_seq);
    disp(['Dual solution Inf relative error (): ' num2str(yopt_rel_err)]);

    disp(['Fraction of optimization solutions: ' ...
        num2str(out.n_opt_solution / length(out.t_opt_seq))]);
    disp(['Fraction of nonempty J_2: ' ...
        num2str(out.n_J2e / length(out.t_opt_seq))]);
end

