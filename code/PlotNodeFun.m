function PlotNodeFun(x0_EZ_range1, BeAwayStabilityMatrix, PushAwayStabilityMatrix, node_index, save_path)

    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end

    % extracting data from the current data
    y = squeeze(BeAwayStabilityMatrix(:, :, node_index));
    z = squeeze(PushAwayStabilityMatrix(:, :, node_index));

    % generate x-axis data for each parameter combination
    x_y = x0_EZ_range1;
    x_z = x0_EZ_range1;

    % create graphs
    figure;
    plot(x_y, y, '-o', 'MarkerSize', 5, 'LineWidth', 1.5); %draw a blue line
    hold on;
    plot(x_z, z, '--', 'Color', 'r', 'LineWidth', 2); %draw a red line

    % add legend
    legend('Node be pushed away from stable state', 'Node pushes network from stable state', 'Location', 'best');

    % add captions and labels
    title(['Node ' num2str(node_index) ' Function'], 'FontSize', 14);
    xlabel('x0(EZ)', 'FontSize', 12);
    ylabel('Node Function', 'FontSize', 12);
    grid on; 
    set(gca, 'FontSize', 12); 

    % adjust the range of axes
    xlim([min(x0_EZ_range1) max(x0_EZ_range1)]);
    ylim([min([y(:); z(:)]) - 10, max([y(:); z(:)]) + 10]); 

    hold off;

    % save figure
    saveas(gcf, fullfile(save_path, ['Node_' num2str(node_index) '.jpg']));
    close(gcf); 
end
