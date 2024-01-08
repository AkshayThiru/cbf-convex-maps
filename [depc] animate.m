function [] = animate(arr, t, x, dt)
    % Animate robots.
    figure;
    hold on
    lims = 15.0;
    xlim([-lims lims]);
    ylim([-lims lims]);
    axis('equal');
    pts = cell(1, length(arr));
    for i = 1:length(arr)
        points_i = arr{i}.geo.plot_outline(x{i}(:,1));
        pts{i} = plot(points_i(1,:), points_i(2,:));
    end
    for k = 1:length(t)
        for i = 1:length(arr)
            points_i = arr{i}.geo.plot_outline(x{i}(:,k));
            pts{i}.XData = points_i(1,:);
            pts{i}.YData = points_i(2,:);
        end
        pause(dt);
    end
end