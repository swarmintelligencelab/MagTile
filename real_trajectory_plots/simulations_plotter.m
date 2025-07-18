clc;
clear;
close all;


data = readtable('Trajectory_L_1.csv');  % Your data
data = data(3000:6000,:);

x = data{:,1};
y = data{:,2};

% Grid parameters
x_min = -30; x_max = 30;
y_min = -30; y_max = 30;
grid_size = 15;

% Interior indices: 3 through 13 (excluding outer 2-cell boundary)
interior_min = 3;
interior_max = grid_size - 2;  % 13

% Compute cell size for interior grid (11 steps)
cell_w = (x_max - x_min) / (grid_size - 4);  % 11 columns
cell_h = (y_max - y_min) / (grid_size - 4);  % 11 rows

% Map x, y to grid indices from 3 to 13
col = floor((x - x_min) ./ cell_w) + 3;
row = floor((y - y_min) ./ cell_h) + 3;

% Clamp to interior region [3, 13]
col = min(interior_max, max(interior_min, col));
row = min(interior_max, max(interior_min, row));

% Remove consecutive duplicates
trajectory = [row, col];
trajectory = trajectory([true; any(diff(trajectory), 2)], :);

% Initialize output path
resolved = [];

% Brute-force walk: one step at a time, strictly row or col
i = 1;
start_pos = trajectory(1,:);
while i < size(trajectory,1)
    end_pos = trajectory(i+1, :);

    curr = start_pos;
    resolved = [resolved; curr];  % add starting point

    while norm([curr(1) - end_pos(1), curr(2) - end_pos(2)]) > 1
        step = [0, 0];

        % Move in row direction first
        if curr(1) ~= end_pos(1)
            step(1) = sign(end_pos(1) - curr(1));
        elseif curr(2) ~= end_pos(2)
            step(2) = sign(end_pos(2) - curr(2));
        end

        curr = curr + step;
        resolved = [resolved; curr];
    end
    start_pos = curr;
    i = i+1;
end

% Add the final point if not already added
if ~isequal(resolved(end,:), trajectory(end,:))
    resolved = [resolved; trajectory(end,:)];
end

% Verify no diagonal or multi-cell jumps remain
steps = diff(resolved);
invalid = any(abs(steps) > 1, 2) | all(abs(steps) == 1, 2);
if any(invalid)
    error('Diagonal or multi-cell jumps remain after resolution!');
end

% Save result (optional: subtract 1 to make 0-indexed)
csvwrite('trajectory_rows_cols.csv', [resolved(:,1)-1, resolved(:,2)-1]);


% ------------------- Plot Real vs Mapped Trajectory -------------------

% Convert grid cell centers to physical coordinates for resolved trajectory
x_resolved = x_min + (resolved(:,2) - 3 + 0.5) * cell_w; % col -> x
y_resolved = y_min + (resolved(:,1) - 3 + 0.5) * cell_h; % row -> y

% Setup background grid points
[xg, yg] = meshgrid(linspace(x_min, x_max, grid_size), linspace(y_min, y_max, grid_size));

figure;
hold on;
axis equal;

% Plot grid as small gray dots
scatter(xg(:), yg(:), 20, 'filled', ...
    'MarkerFaceColor', [0.6 0.6 0.6], ...
    'HandleVisibility', 'off');
% Plot resolved (mapped) trajectory in light orange with squares
plot(x_resolved, y_resolved, '-s', ...
    'Color', [1, 0.6, 0], ...
    'MarkerFaceColor', [1, 0.6, 0], ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Mapped Trajectory');

% Plot real (continuous) trajectory in shaded light blue
plot(x, y, '-', ...
    'Color', [0.3 0.75 0.93], ...
    'LineWidth', 2, ...
    'DisplayName', 'Real Trajectory');

% Mark start point (real)
scatter(x(1), y(1), 120, 'o', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [1, 0.6, 0], ...
    'DisplayName', 'Start');

% Mark end point (real) in purple
scatter(x(end), y(end), 120, 'o', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.5 0 0.5], ...
    'DisplayName', 'End');
% Axis formatting
x_ticks = linspace(x_min, x_max, 5);
y_ticks = linspace(y_min, y_max, 5);
xlim([x_min, x_max]);
ylim([y_min, y_max]);
xlabel('$\mathcal{X}$ (cm)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\mathcal{Y}$ (cm)', 'Interpreter', 'latex', 'FontSize', 16);

% Custom tick labels (optional for matching figure)
xticks(linspace(x_min, x_max, 5));
yticks(linspace(y_min, y_max, 5));
text(x_ticks(1), y_min - 2, '$-x_{\max}$', 'Interpreter','latex', 'HorizontalAlignment','center', 'FontSize', 14);
text(x_ticks(3), y_min - 2, '0',           'Interpreter','latex', 'HorizontalAlignment','center', 'FontSize', 14);
text(x_ticks(5), y_min - 2, '$x_{\max}$',  'Interpreter','latex', 'HorizontalAlignment','center', 'FontSize', 14);

% Y-tick LaTeX labels
text(x_min - 3, y_ticks(1), '$-y_{\max}$', 'Interpreter','latex', 'HorizontalAlignment','right', 'FontSize', 14);
text(x_min - 3, y_ticks(3), '0',           'Interpreter','latex', 'HorizontalAlignment','right', 'FontSize', 14);
text(x_min - 3, y_ticks(5), '$y_{\max}$',  'Interpreter','latex', 'HorizontalAlignment','right', 'FontSize', 14);

% Legend
legend('Location', 'northeastoutside', 'Interpreter', 'latex');

box on;