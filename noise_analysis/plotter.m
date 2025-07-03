clear all
close all
clc

data = load("data_100.mat");

% Example: Replace this with your real signal matrix
% E should be [n_trials x n_timepoints]
n_trials = 100;
p = data.prob_;  % Time from 0 to 2 seconds
counter = data.counter_';

% Compute statistics
E_mean = mean(counter, 1);  % Mean across trials
E_std  = std(counter, 0, 1);  % Standard deviation across trials

% Plot
% figure; hold on; box on;
% Create figure with transparent background
figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 4 3]); 
axes('Color', 'w'); hold on;

% Plot shaded standard deviation
fill([p, fliplr(p)], ...
     [E_mean + E_std, fliplr(E_mean - E_std)], ...
     [0.6 0.8 1], ...              % light blue color
     'FaceAlpha', 0.5, ...
     'EdgeColor', 'none');

% Plot mean line in bold blue
plot(p, E_mean, 'b', 'LineWidth', 2);

% Labels
xlabel('$p$', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'k');
ylabel('$N_{\mathrm{A}}$', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'k');

% Legend
legend({'$\bar{N_{\mathrm{A}}} \pm \mathrm{SD}$', '$\bar{N_{\mathrm{A}}}$'}, ...
       'Interpreter', 'latex', 'Location', 'southeast', 'Color', 'w', 'TextColor', 'k');

% Axes formatting
set(gca, 'FontSize', 12);
xlim([0, 1]);
ylim padded;

set(gca, 'LooseInset', [0 0 0 0]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca, ...
    'XColor', 'k', ...  % X-axis tick color
    'YColor', 'k', ...  % Y-axis tick color
    'TickLabelInterpreter', 'latex', ...
    'FontSize', 12);




% Example: Replace this with your real signal matrix
% E should be [n_trials x n_timepoints]
n_trials = 100;
p = data.prob_;  % Time from 0 to 2 seconds
e_ = data.e_;
e_ = reshape(e_, 100, 20).';
e_ = e_';

% Compute statistics
E_mean = mean(e_, 1);  % Mean across trials
E_std  = std(e_, 0, 1);  % Standard deviation across trials

% Plot
% figure; hold on; box on;
% Create figure with transparent background
figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 4 3]); 
axes('Color', 'w'); hold on;

% Plot shaded standard deviation
fill([p, fliplr(p)], ...
     [E_mean + E_std, fliplr(E_mean - E_std)], ...
     [0.6 0.8 1], ...              % light blue color
     'FaceAlpha', 0.5, ...
     'EdgeColor', 'none');

% Plot mean line in bold blue
plot(p, E_mean, 'b', 'LineWidth', 2);

% Labels
xlabel('$p$', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'k');
ylabel('$E$', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'k');

% Legend
legend({'$\bar{E} \pm \mathrm{SD}$', '$\bar{E}$'}, ...
       'Interpreter', 'latex', 'Location', 'southeast', 'Color', 'w', 'TextColor', 'k');

% Axes formatting
set(gca, 'FontSize', 12);
xlim([0, 1]);
ylim padded;

set(gca, 'LooseInset', [0 0 0 0]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca, ...
    'XColor', 'k', ...  % X-axis tick color
    'YColor', 'k', ...  % Y-axis tick color
    'TickLabelInterpreter', 'latex', ...
    'FontSize', 12);
