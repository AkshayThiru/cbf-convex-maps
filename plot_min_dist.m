function [] = plot_min_dist(T, t, arr, log, const, save_figures)
    % Set colours.
    color_band = [[0, 0.4470, 0.7410]
        [0.8500, 0.3250, 0.0980];
        [0.9290, 0.6940, 0.1250];
        [0.4940, 0.1840, 0.5560];
        [0.4660, 0.6740, 0.1880];
        [0.6350, 0.0780, 0.1840]];
    
    pair_idx = @(i,j) ((2*length(arr)-i)*(i-1))/2 + j-i; % j > i.
    
    fig_L = 600;
    fig_H = 200;
    
    %% Plot minimum distance.
    figure('Renderer', 'painters', 'Position', [0 0 4/3*fig_L fig_H]);
    til = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
    
    % Minimum distance between robot_idx and other robots.
    robot_idx = 1;
    nexttile;
    hold on;
    for i = setdiff(1:length(arr), robot_idx)
        k = pair_idx(min(i, robot_idx), max(i, robot_idx));
        plot(t(1:end-1), log.h(:,k), '-', 'Color', color_band(i,:), 'LineWidth', 1.5);
    end
    plot([t(1), t(end)], const.e*[1, 1], '--', 'Color', 'k', 'LineWidth', 1.0);
    ax = gca;
    ax.XTick = 0:5:T;
    ax.YTick = [10^-3, 10^-2, 10^-1, 1, 10, 100];
    xlim([0 T]);
    ylim([10^-2 300]);
    box on
    set(gca,'LineWidth', 1, 'FontSize', 10);
    set(gca, 'YScale', 'log');
    legend_text = cell(length(arr),1);
    for i = 1:length(arr)-1
        idx_set = setdiff(1:length(arr), robot_idx);
        legend_text{i} = strcat('$h^{',num2str(robot_idx),num2str(idx_set(i)),'}(t)$');
    end
    legend_text{length(arr)} = '$\epsilon_1^2$';
    h_l = flipud(get(gca,'Children'));
    h_legend = legend(h_l, legend_text);
    h_legend.FontSize = 12;
    h_legend.ItemTokenSize = [15, 18];
    h_legend.NumColumns = 1;
    h_legend.Location = 'southeast';
    set(h_legend, 'Interpreter','latex');
    hold off;
    
    % Plot min(h^{ij}).
    nexttile;
    hold on;
    plot(t(1:end-1), min(log.h, [], 2), '-', 'Color', color_band(end,:), 'LineWidth', 1.5);
    plot([t(1), t(end)], const.e*[1, 1], '--', 'Color', 'k', 'LineWidth', 1.0);
    ax = gca;
    ax.XTick = 0:5:T;
    ax.YTick = [10^-3, 10^-2, 10^-1, 1, 10, 100];
    set(gca, 'YTickLabel', []);
    xlim([0 T]);
    ylim([10^-2 300]);
    box on
    set(gca,'LineWidth', 1, 'FontSize', 10);
    set(gca, 'YScale', 'log');
    h_l = flipud(get(gca,'Children'));
    h_legend = legend(h_l, {'$\min_{i,j}(h^{ij})$', '$\epsilon_1^2$'});
    h_legend.FontSize = 12;
    h_legend.ItemTokenSize = [15, 18];
    h_legend.NumColumns = 1;
    h_legend.Location = 'northeast';
    set(h_legend, 'Interpreter','latex');
    hold off;
    
    xlabel(til, 'Time (s)','interpreter','latex','FontSize', 15);
    ylabel(til, '(m$^2$)','interpreter','latex','FontSize', 15);
    
    if save_figures
        print(gcf, strcat('./fig/strict-convex-h.png'), '-dpng', '-r1000');
        print(gcf, strcat('./fig/strict-convex-h.eps'), '-depsc', '-r600');
    end
end
     