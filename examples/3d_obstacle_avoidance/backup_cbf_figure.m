function [] = backup_cbf_figure()
    %% Backup trajectory + backup CBF figure.
    % Set figure.
    % Manually set the size to (250, 150) points.
    fig = figure();
    % Create tiledlayout.
    til = tiledlayout(fig, 1, 1); % (rows, columns).
    til.TileSpacing = "tight";
    til.Padding = "tight";
    ax = nexttile;

    % Set polytopes.
    V1 = [0, 1, 0; 1, 0, 0; -1, 0, 0];
    V2 = [1, 1, 0; 1, -1, 0; -1, -1, 0; -1, 1, 0] + ones(4, 1) * [3, 3, 0];
    floor = [-2, -1, -0.01; 7, -1, -0.01; 7, 5, -0.01; -2, 5, -0.01];
    hold on
    fill3(ax, V1(:, 1), V1(:, 2), V1(:, 3), [200, 0, 0]/255, ...
        V2(:, 1), V2(:, 2), V2(:, 3), [0, 200, 0]/255, ...
        'FaceAlpha', 1, 'EdgeColor', 'k');
    fill3(floor(:, 1), floor(:, 2), floor(:, 3), 'k', ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none');
    % Set displacements.
    disp1 = [1, 2, 4];
    disp2 = [2, -2, 4];
    patch1 = extrude_2d_patch(V1, disp1);
    patch2 = extrude_2d_patch(V2, disp2);
    for i = 1:length(patch2)
        Vi = patch2{i};
        fill3(ax, Vi(:, 1), Vi(:, 2), Vi(:, 3), 'g', ...
            'FaceAlpha', 0.1, 'EdgeColor', 'k');
    end
    for i = 1:length(patch1)
        Vi = patch1{i};
        fill3(ax, Vi(:, 1), Vi(:, 2), Vi(:, 3), 'r', ...
            'FaceAlpha', 0.1, 'EdgeColor', 'k');
    end
    % Draw arrows.
    quiver3(0, 0.5, 0, disp1(1), disp1(2), 0, 0.75, ...
        'Color', 'k', 'LineWidth', 1.25, 'MaxHeadSize', 1);
    quiver3(3, 3, 0, disp2(1), disp2(2), 0, 0.75, ...
        'Color', 'k', 'LineWidth', 1.25, 'MaxHeadSize', 1);
    % Compute minimum distance.
    P1 = StaticPolytope([V1; V1 + ones(size(V1, 1), 1) * disp1]);
    P2 = StaticPolytope([V2; V2 + ones(size(V2, 1), 1) * disp2]);
    [~, z_opt, ~,~,~] = minimum_distance(ones(0, 1), P1, ones(0, 1), P2);
    % plot3(z_opt(1), z_opt(2), z_opt(3), 'ok');
    % plot3(z_opt(4), z_opt(5), z_opt(6), 'ok');
    % Add text.
    text(0, 0, 6, '$\mathcal{P}^i$', 'Interpreter', 'latex', ...
        'FontUnits', 'points', 'FontSize', 10);
    text(0, 0, 6, '$\bar{\mathcal{P}}^i$', 'Interpreter', 'latex', ...
        'FontUnits', 'points', 'FontSize', 10);
    text(0, 0, 6, '$\mathcal{P}^j$', 'Interpreter', 'latex', ...
        'FontUnits', 'points', 'FontSize', 10);
    text(0, 0, 6, '$\bar{\mathcal{P}}^j$', 'Interpreter', 'latex', ...
        'FontUnits', 'points', 'FontSize', 10);
    text(0, 0, 6, '$z$', 'Interpreter', 'latex', ...
        'FontUnits', 'points', 'FontSize', 10);
    text(0, 0, 6, '$\xi$', 'Interpreter', 'latex', ...
        'FontUnits', 'points', 'FontSize', 10);
    % Axis properties.
    % ax.XLim = [-2, 7];
    % ax.YLim = [-1, 5];
    % ax.ZLim = [-0.01, 2];
    ax.XTick = []; ax.YTick = []; ax.ZTick = [];
    ax.Box = "off";
    view([24, 26]);
    axis equal

    %% Softmax over-approximation.
    fig = figure();
    % Create tiledlayout.
    til = tiledlayout(fig, 1, 2); % (rows, columns).
    til.TileSpacing = "tight";
    til.Padding = "tight";
    til.TileIndexing = "rowmajor"; % (default).

    % Set polytope and smooth sets.
    P = StaticPolytope([0, 1, 0; 1, -1, 0; -1, -1, 0; 0, 0, 1]);
    R_ = eye(3);
    x = [0; 0; 0; R_(:)];
    C1 = SmoothPolytope(3, P, 1);
    C1.center = [0; 0; 0]; C1. radius = 10;
    C2 = SmoothPolytope(3, P, 5);
    C2.center = [0; 0; 0]; C2. radius = 10;

    % Plot.
    ax = nexttile;
    P.plot_surf([], ax, 'c', 1, []);
    C1.plot_surf(x, ax, 'r', 0.25, 0.5);
    ax.XTick = []; ax.YTick = []; ax.ZTick = [];
    ax.Box = "off";
    axis equal
    % Set legends.
    legends = {"$\kappa = 1$"};
    child_ = flipud(ax.Children);
    leg = legend([child_(2)], legends); % Each surface has 10 plots.
    leg.FontSize = 10;
    leg.ItemTokenSize = [25, 25]; % Box size, default [30, 18], *Undocumented*.
    leg.NumColumns = 1;
    leg.Location = "northeast";
    leg.Interpreter = "latex";
    % legend off;
    leg.BoxFace.ColorType = 'truecoloralpha';
    leg.BoxFace.ColorData = uint8(255 * [1, 1, 1, 0.5]');

    ax = nexttile;
    P.plot_surf([], ax, 'c', 1, []);
    C2.plot_surf(x, ax, 'r', 0.25, 0.5);
    ax.XTick = []; ax.YTick = []; ax.ZTick = [];
    ax.Box = "off";
    axis equal
    % Set legends.
    legends = {"$\kappa = 5$"};
    child_ = flipud(ax.Children);
    leg = legend([child_(2)], legends); % Each surface has 10 plots.
    leg.FontSize = 10;
    leg.ItemTokenSize = [15, 15]; % Box size, default [30, 18], *Undocumented*.
    leg.NumColumns = 1;
    leg.Location = "northeast";
    leg.Interpreter = "latex";
    % legend off;
    leg.BoxFace.ColorType = 'truecoloralpha';
    leg.BoxFace.ColorData = uint8(255 * [1, 1, 1, 0.5]');
end

function [patches] = extrude_2d_patch(V, disp)
    % Extrude V along disp.
    % Returns E+2 patches, where E is the number of edges of the 2d patch V.
    E = size(V, 1);
    patches = cell(E + 2, 1);
    patches{1} = V;
    patches{E+2} = V + ones(E, 1) * disp;

    idx = [1:E, 1];
    for i = 1:E
        patches{i+1} = [V(idx(i), :);
            V(idx(i+1), :);
            V(idx(i+1), :) + disp;
            V(idx(i), :) + disp];
    end
end
