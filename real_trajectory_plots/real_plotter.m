clc;
clear;
close all;

% Load reference (grid-based) trajectory [row, col]
data_2 = readtable('Trajectory_2.csv');
data_2(:, 3) = num2cell((0:0.2:(height(data_2)-1)*0.2)');  % Add time column
row = data_2{:,1};
col = data_2{:,2};
t = data_2{:,3};

% Load real trajectory (e.g., camera/tracker output)
data = readtable('Trajectory_L_1.csv');  
data = data(3000:6000,:);  % Focus on central portion
x_real = data{:,1};
y_real = data{:,2};
t_real = linspace(0, 0.2*(length(x_real)-1), length(x_real));

target_min = -15;
target_max = 15;

% Normalize and rescale
x_real = (x_real - min(x_real)) / (max(x_real) - min(x_real));  % [0,1]
x_real = x_real * (target_max - target_min) + target_min;       % [-15,15]

y_real = (y_real - min(y_real)) / (max(y_real) - min(y_real));  % [0,1]
y_real = y_real * (target_max - target_min) + target_min;       % [-15,15]

% Downsample real trajectory to match reference
n_ref = height(data_2);
x_real_ds = x_real(round(linspace(1, length(x_real), n_ref)));
y_real_ds = y_real(round(linspace(1, length(y_real), n_ref)));
t_real_ds = t;  % Use same time vector for simplicity

% Grid parameters (must match mapping logic)
x_min = -15; x_max = 15;
y_min = -15; y_max = 15;
grid_size = 15;
interior_offset = 3;
cell_w = (x_max - x_min) / (grid_size - 4);
cell_h = (y_max - y_min) / (grid_size - 4);

% Convert grid indices to physical coordinates
x_ref = x_min + (col - interior_offset + 0.5) * cell_w;
y_ref = y_min + (row - interior_offset + 0.5) * cell_h;

% --- 3D Plot ---
fig = figure( ...
    'Color', 'w', ...
    'Units', 'inches', ...
    'Position', [1 1 4 3], ...
    'Toolbar', 'figure', ...
    'MenuBar', 'figure' ...
);
ax = axes('Parent', fig, 'Color', 'w'); 
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on'); axis(ax, 'equal');

% Real trajectory: solid black
h1 = plot3(t_real_ds, x_real_ds, y_real_ds, 'r--', 'LineWidth', 2);

% Reference trajectory: dashed red
h2 = plot3(t, x_ref, y_ref, 'k-', 'LineWidth', 2);

yticks([-10, 10])
zticks([-10, 10])

% yticks([-10,10])
% zticks([-10,10])
% Labels
xlabel('$t$ (s)', 'Interpreter','latex', 'FontSize',14);
ylabel('$\mathcal{X}$ (cm)', 'Interpreter','latex', 'FontSize',14);
zlabel('$\mathcal{Y}$ (cm)', 'Interpreter','latex', 'FontSize',14);

% Set correct 3D view to match your sample figure
view([135, 30]);  % [azimuth, elevation]

% Axis formatting
set(gca, ...
    'FontSize', 12, ...
    'TickLabelInterpreter', 'latex', ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'ZColor', 'k');




