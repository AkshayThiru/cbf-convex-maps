clear

current_path = path;
path(path, genpath('..\'));


%% Define convex sets.
V = [-1, -1; -1, 1; 1, 1; 1, -1; 0, 0; 1, 1; 0, 2];
V = [V, -ones(7, 1); V, ones(7, 1)];
Q = [3, 1, 0; 1, 2, 0; 0, 0, 1];
sp = StaticPolytope(V);
smp = SmoothPolytope(3, sp, 10);
re = RigidEllipsoid(3, Q);

x_sp = zeros(0, 1);
p_smp = [10; 10; 10];
theta_z = 15 * pi/180;
R_smp = [cos(theta_z), -sin(theta_z), 0;
    sin(theta_z), cos(theta_z), 0;
    0, 0, 1];
x_smp = [p_smp; reshape(R_smp, [9, 1])];
p_re = [0; 0; 0];
R_re = eye(3);
x_re = [p_re; reshape(R_re, [9, 1])];


%% Minimum distance test.

% Distance between SmoothPolytope and StaticPolytope.
[dist2, z_opt, y_opt, ~, ~] = minimum_distance(x_sp, sp, x_smp, smp, ...
    [], 'interior-point');
solution_optimality(x_sp, sp, x_smp, smp, z_opt, y_opt, dist2);
% Run time comparision between 'interior-point' and 'sqp'.
[~,~,~,~, out_ipopt] = minimum_distance(x_sp, sp, x_smp, smp, ...
    [], 'interior-point');
[~,~,~,~, out_sqp] = minimum_distance(x_sp, sp, x_smp, smp, ...
    [], 'sqp');
disp(['IPOpt run time (s): ' num2str(out_ipopt.run_time)]);
disp(['SQP run time (s)  : ' num2str(out_sqp.run_time)]);
disp(' ');

% Distance between SmoothPolytope and RigidEllipsoid.
[dist2, z_opt, y_opt, ~, ~] = minimum_distance(x_smp, smp, x_re, re, ...
    [], 'interior-point');
solution_optimality(x_smp, smp, x_re, re, z_opt, y_opt, dist2);
% Run time comparision between 'interior-point' and 'sqp'.
[~,~,~,~, out_ipopt] = minimum_distance(x_smp, smp, x_re, re, ...
    [], 'interior-point');
[~,~,~,~, out_sqp] = minimum_distance(x_smp, smp, x_re, re, ...
    [], 'sqp');
disp(['IPOpt run time (s): ' num2str(out_ipopt.run_time)]);
disp(['SQP run time (s)  : ' num2str(out_sqp.run_time)]);
disp(' ');


%%
path(current_path);


%% Function for checking solution optimality.
function [] = solution_optimality(x1, C1, x2, C2, z_opt, y_opt, dist2)
    z1 = z_opt(1:C1.nz); z2 = z_opt(C1.nz+1:end);
    y1 = y_opt(1:C1.nr); y2 = y_opt(C1.nr+1:end);
    [A1, ~, dAdz1, ~,~] = C1.derivatives(x1, z1, zeros(C1.nr, 1));
    [A2, ~, dAdz2, ~,~] = C2.derivatives(x2, z2, zeros(C2.nr, 1));
    % Stationarity check:
    stationarity = [2 * (z1 - z2) + dAdz1' * y1;
        2 * (z2 - z1) + dAdz2' * y2];
    disp(['Stationarity error: ' num2str(norm(stationarity, Inf))]);
    % Objective value check:
    obj_value = (z1 - z2)' * (z1 - z2) - dist2;
    disp(['Objective value error: ' num2str(obj_value)]);
    % Primal feasibility check:
    primal_feasibility = [max(A1); max(A2)];
    disp(['Primal feasibility: ' num2str(primal_feasibility')]);
    % Dual feasibility check:
    dual_feasibility = [min(y1); min(y2)];
    disp(['Dual feasibility: ' num2str(dual_feasibility')]);
    % Complementary slackness check:
    complementary_slackness = [y1' * A1; y2' * A2];
    disp(['Complementary slackness: ' num2str(complementary_slackness')]);
end
