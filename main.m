path(path, genpath('.\'));

global safety
safety = true;

%% Initialize robots.
arr = init_robots();

%% Set constants.
const.N = (length(arr)*(length(arr)-1))/2;
const.m = ones(length(arr)+1, 1);
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
statistics(time, log, const);

% Plot trajectories.
figure(1)
hold on
for i = 1:length(arr)
    plot(x{i}(1,:), x{i}(2,:));
end

% Plot states.
figure(2)
hold on
arr_idx = 2; state_idx = 3;
plot(t, x{arr_idx}(state_idx,:));

% Plot minimum distance.
figure(3)
hold on
arr_idx = 4;
plot(t(1:end-1), log.h(:,arr_idx));

% Animate robots.
animate(arr, t, x, dt);

%% Error handling.
if ~isempty(ME)
   rethrow(ME);
end
