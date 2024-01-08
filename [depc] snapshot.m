function [] = snapshot(ts, x, arr, save_snapshot)
    % Set colours.
    color_band = [[0, 0.4470, 0.7410]
        [0.8500, 0.3250, 0.0980];
        [0.9290, 0.6940, 0.1250];
        [0.4940, 0.1840, 0.5560];
        [0.4660, 0.6740, 0.1880];
        [0.6350, 0.0780, 0.1840]];
    
    figure('Renderer', 'painters', 'Position', [0 0 600 300]);
    set(gca,'LooseInset',get(gca,'TightInset'));
    hold on;
    
    % Plot robots.
    for i = 1:length(arr)
        for k = 1:length(ts)
            points_i = arr{i}.geo.plot_outline(x{i}(:,ts(k)));
            % Fill interior.
            fill(points_i(1,:), points_i(2,:), color_band(i,:), ...
                'FaceAlpha', 1 - (k-0)/(2*length(ts)-0), 'EdgeColor', 'k');
        end
        % Plot robot trajectories.
        plot(x{i}(1,:), x{i}(2,:), '--', 'Color', color_band(i,:), 'LineWidth', 1.0);
        % Plot initial and final positions.
        plot(arr{i}.xi(1), arr{i}.xi(2), 'ok');
        plot(arr{i}.xf(1), arr{i}.xf(2), '*k');
    end
    xlim([-11, 11]);
    ylim([-5.5, 5.5]);
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    title('Snapshots of obstacle avoidance for strictly convex sets', ...
        'interpreter', 'latex', 'FontSize', 15);
    axis('equal');
    box off;
    axis off;
    grid off;
    hold off;
    
    if save_snapshot
        print(gcf, './fig/strict-convex-snapshot.png', '-dpng', '-r1000');
        print(gcf, './fig/strict-convex-snapshot.eps', '-depsc', '-r600');
    end
end
