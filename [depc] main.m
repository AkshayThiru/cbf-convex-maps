path(path, genpath('.\'));

global safety
safety = true;

%% Display options.
display_solve_stats = true;
plot_stat_hist      = false;
plot_figures        = true;
save_figures        = false;
plot_snapshot       = true;
save_snapshot       = false;
display_animation   = true;

%% Initialize robots.
arr = init_robots();

%% Set constants.
const.a_cbf = 1; % ECBF rate.
const.e     = 0.1; % Safety margin [m^2].
const.M_L   = 1e3; % Bound on dL. |dL|.|dA| ~ |dh/dt|.

const.N = (length(arr)*(length(arr)-1))/2; % NC2.
const.m = ones(length(arr)+1, 1); % Index start and stop for inputs.
for i = 1:length(arr)
    const.m(i+1) = const.m(i) + arr{i}.m;
end
const.M = const.m(end) - 1;

%% Initialize data structures.
T = 20;
tspan = [0 T];
dt = 0.05;

t = (0:dt:T+dt)';
x = cell(1, length(arr));
for i = 1:length(x)
    x{i} = zeros(arr{i}.n, size(t,1));
    x{i}(:,1) = arr{i}.x;
end
time = zeros(size(t,1),1);
prev = cell(const.N, 1);
for i = 1:length(arr)
    for j = i+1:length(arr)
        k = ((2*length(arr)-i)*(i-1))/2 + j-i;
        prev{k} = zeros(2*arr{1}.l, 1);
    end
end
log = struct('h', [], 't_opt', []);
debug = [];
display_text = 0;
ME = [];

%% Run simulation.
for k = 1:length(t)-1
    % Evaluate control input.
    try
        start = tic;
        [u_t, prev, log_t, debug_t] = control(arr, prev, const);
        time(k) = toc(start);
    catch ME
        t = t(1:k);
        for i = 1:length(arr)
            x{i} = x{i}(:,1:k);
        end
        time = time(1:k);
        break;
    end
    % Statistics.
    fprintf(repmat('\b', 1, display_text));
    display_text = fprintf("time: %2.2f s, loop time: %1.4f s, frequency: %3.1f Hz", ...
        t(k), time(k), 1/time(k));
    % Simulate system.
    for i = 1:length(arr)
        [~,x_te] = ode45(@(t,x) arr{i}.dyn(x, u_t{i}), [0 dt], x{i}(:,k));
        x{i}(:,k+1) = x_te(end,:)';
        arr{i}.x = x_te(end,:)';
    end
    % Logging.
    logn = fieldnames(log);
    for i = 1:numel(logn)
        log.(logn{i})(k,:) = log_t.(logn{i});
    end
end
fprintf(repmat('\b', 1, display_text));

% Average statistics.
if display_solve_stats
    statistics(time, log, const, plot_stat_hist);
end

% Plot minimum distance.
if plot_figures
    plot_min_dist(T, t, arr, log, const, save_figures);
end

% Plot snapshot.
if plot_snapshot
    ts = [0 8];
    t_idx = round(ts/dt) + 1;
    snapshot(t_idx, x, arr, save_snapshot);
end

% Animate robots.
if display_animation
    animate(arr, t, x, dt);
end

%% Error handling.
if ~isempty(ME)
   rethrow(ME);
end
