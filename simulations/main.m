clc;
clear;
close all;

% =========================================
% Constants and Parameters
% =========================================

xGridSize = 15;  % Grid dimensions for static dipoles (x-direction)
yGridSize = 15;  % Grid dimensions for static dipoles (y-direction)
GridSpacing = 2.159;  % Spacing between static dipoles: 2.159cm
H = 2.485; % Height of the plane above the origin: 2.485 cm

alpha = 22042.7968;  % Magnetic force coefficient: 22042.7968 calculated in cm
beta = 15.5488; % Damping coefficient: 15.5488 calculated in cm

MagRng = 1.5*GridSpacing; % Magnetic force range
MagForceFar = 1; % Magnetic force multiplier: 1
MagRngTest = 3*MagRng; % Magnetic disk distance test range
MagRngWeight = 2*MagRng; % Weight adjustion range

dt = 0.001;  % Simulation time step
T = 200; % Total simulation time
vt = 0:dt:T;  % Simulation time step vector
d = 10;  % Sampling period
Td = d/dt; % Time steps per sampling period

% =========================================
% Generate Grid for Static Dipoles
% =========================================

[xGrid, yGrid] = meshgrid(linspace(-(xGridSize-1)/2, (xGridSize-1)/2, xGridSize)*GridSpacing, ...
    linspace((yGridSize-1)/2, -(yGridSize-1)/2, yGridSize)*GridSpacing);
zGrid = zeros(size(xGrid));  % Assume all dipoles are at z = 0

% Initialize an mxm cell matrix for coordinates
xyGridCell = cell(xGridSize, yGridSize);
for i = 1:xGridSize
    for j = 1:yGridSize
        xyGridCell{i, j} = [xGrid(i, j), yGrid(i, j)]';
    end
end

% =========================================
% Create Network (Adjacency Matrix)
% =========================================

numCoil = numel(xGrid);

% Define the connectivity in a grid (8 neighbors: 4 nearest + 4 diagonal)
A = zeros(numCoil, numCoil);  % Initialize the adjacency matrix
for i = 1:xGridSize
    for j = 1:yGridSize
        currentInd = sub2ind(size(xGrid), i, j);
        neighbors = [i, j-1; i, j+1; i-1, j; i+1, j; ...
            i-1, j-1; i-1, j+1; i+1, j-1; i+1, j+1];
        for k = 1:length(neighbors)
            ni = neighbors(k, 1);
            nj = neighbors(k, 2);
            if ni >= 1 && ni <= xGridSize && nj >= 1 && nj <= yGridSize
                neighborInd = sub2ind(size(xGrid), ni, nj);
                A(currentInd, neighborInd) = norm([i,j]-[ni, nj]);
                A(neighborInd, currentInd) = norm([i,j]-[ni, nj]);
            end
        end
    end
end

% =========================================
% Disk and Reference Initialization
% =========================================

% Disk 1
x1Disk = zeros(size(vt));
y1Disk = zeros(size(vt));
vx1Disk = zeros(size(vt));
vy1Disk = zeros(size(vt));

% Reference 1 orbits of trajectory and control
xy1Ref0Ind = [34, 35, 51, 66, 80, 79, 63, 48]; % Ctr diamond
xy1Ref0 = [xGrid(xy1Ref0Ind); yGrid(xy1Ref0Ind)]; % In location

u1Ref0Ind = xy1Ref0Ind; % Same as trajectory

% Initial positions and velocities of disk 1 (moving dipole)
x1Disk(1) = 10;
y1Disk(1) = -10;
z1Disk = H;
vx1Disk(1) = 0;  % Initial velocity in x-direction
vy1Disk(1) = 0;  % Initial velocity in y-direction

% Reference 1 signals of trajectory and control
xy1Ref = repmat(xy1Ref0, 1, ceil(T/d));  % Reference signal for trajectory
xy1Ref = xy1Ref(:, 1:T/d+1);
u1RefInd = repmat(u1Ref0Ind, 1, ceil(T/d)); % Reference signal for control
u1RefInd = u1RefInd(:, 1:T/d+1);

% Reference 1 signals in sim steps
xy1RefPlot = [kron(xy1Ref(:, 1:end-1), ones(1, Td)), xy1Ref(:, end)];
u1RefIndPlot = [kron(u1RefInd(:, 1:end-1), ones(1, Td)), u1RefInd(:, end)];

% Disk 2
x2Disk = zeros(size(vt));
y2Disk = zeros(size(vt));
vx2Disk = zeros(size(vt));
vy2Disk = zeros(size(vt));

% Reference 2 orbits of trajectory and control
xy2Ref0Ind = [34, 35, 51, 66, 80, 79, 63, 48] + 8*15;
xy2Ref0 = [xGrid(xy2Ref0Ind); yGrid(xy2Ref0Ind)]; % In location

u2Ref0Ind = xy2Ref0Ind; % Same as trajectory

% Initial positions and velocities of disk 2 (moving dipole)
x2Disk(1) = -8;
y2Disk(1) = -10;
z2Disk = H;
vx2Disk(1) = 0;  % Initial velocity in x-direction
vy2Disk(1) = 0;  % Initial velocity in y-direction

% Reference 2 signals of trajectory and control
xy2Ref = repmat(xy2Ref0, 1, ceil(T/d));  % Reference signal for trajectory
xy2Ref = xy2Ref(:, 1:T/d+1);
u2RefInd = repmat(u2Ref0Ind, 1, ceil(T/d)); % Reference signal for control
u2RefInd = u2RefInd(:, 1:T/d+1);

% Reference 1 signals in sim steps
xy2RefPlot = [kron(xy2Ref(:, 1:end-1), ones(1, Td)), xy2Ref(:, end)];
u2RefIndPlot = [kron(u2RefInd(:, 1:end-1), ones(1, Td)), u2RefInd(:, end)];


% =========================================
% Simulation Loop
% =========================================

figure(99);
axis([[-1, 1]*(floor(GridSpacing*(xGridSize-1)/2)+1), ...
    [-1, 1]*(floor(GridSpacing*(yGridSize-1)/2)+1)]);
pbaspect([1, 1, 1]); daspect([1, 1, 1]);
box on; grid on;
hold on;
plot(xGrid(:), yGrid(:), 'g*');

u1IndPlot = zeros(size(vt)); % Control signal in sim steps
e1NormPlot = zeros(size(vt)); % Error signal in sim steps
F1Plot = zeros(2, length(vt)); % Force in sim steps

u2IndPlot = zeros(size(vt)); % Control signal in sim steps
e2NormPlot = zeros(size(vt)); % Error signal in sim steps
F2Plot = zeros(2, length(vt)); % Force in sim steps

for s = 1:T/d % Sampling and control loop

    A2New = A;
    A1New = A;

    disp(['==== Loop ', num2str(s), ' ====']);

    SampleTimeStart = (s-1)*Td + 1;
    SampleTimeEnd = s*Td;
    SampleTimeHalf = (s-1/2)*Td;

    % Observe disk location
    xy1Sample = [x1Disk(SampleTimeStart); y1Disk(SampleTimeStart)];

    % Input reference signals
    xy1RefSample = xy1Ref(:, s);
    u1RefIndSample = u1RefInd(s);

    xy2Sample = [x2Disk(SampleTimeStart); y2Disk(SampleTimeStart)];
    xy2RefSample = xy2Ref(:, s);
    u2RefIndSample = u2RefInd(s);

    % Compute error between disk and reference
    e1Sample = norm(xy1Sample - xy1RefSample);
    e2Sample = norm(xy2Sample - xy2RefSample);
    e12Sample = norm(xy1Sample - xy2Sample);

    MagForce1 = 1;
    MagForce2 = 1;
    if e12Sample > MagRngTest*0.75
        
        if (e1Sample<MagRng) && (e2Sample<MagRng)
            u1IndSample1 = u1RefIndSample;
            u1IndSample2 = u1RefIndSample;
            

            u2IndSample1 = u2RefIndSample;
            u2IndSample2 = u2RefIndSample;
            

        elseif (e1Sample>MagRng) && (e2Sample<MagRng)

            [u1IndSample1, u1IndSample2] = SelectCoilFar(xy1Sample, u1RefIndSample, A1New, ...
                xyGridCell, numCoil, MagRng, s, e1Sample, xy1RefSample, 1);
            % [Debug] Compare using the original weights
            [~, ~] = SelectCoilFar(xy1Sample, u1RefIndSample, A, ...
                xyGridCell, numCoil, MagRng, s, e1Sample, xy1RefSample, 1);
            MagForce1 = MagForceFar;

            u2IndSample1 = u2RefIndSample;
            u2IndSample2 = u2RefIndSample;
            MagForce2 = 1;

        elseif (e1Sample<MagRng) && (e2Sample>MagRng)

            [u2IndSample1, u2IndSample2] = SelectCoilFar(xy2Sample, u2RefIndSample, A2New, ...
                xyGridCell, numCoil, MagRng, s, e2Sample, xy2RefSample, 1);
            % [Debug] Compare using the original weights

            u1IndSample1 = u1RefIndSample;
            u1IndSample2 = u1RefIndSample;
            MagForce1 = 1;

        else

            [u1IndSample1, u1IndSample2] = SelectCoilFar(xy1Sample, u1RefIndSample, A1New, ...
                xyGridCell, numCoil, MagRng, s, e1Sample, xy1RefSample, 1);

            [u2IndSample1, u2IndSample2] = SelectCoilFar(xy2Sample, u2RefIndSample, A2New, ...
                xyGridCell, numCoil, MagRng, s, e2Sample, xy2RefSample, 1);

        end

        u2IndSample1_prev = u2IndSample1;
        u2IndSample2_prev = u2IndSample2;

        u1Sample1 = zeros(numCoil, 1); % Coil vector
        u1Sample2 = zeros(numCoil, 1);
        u1Sample1(u1IndSample1) = MagForce1;
        u1Sample2(u1IndSample2) = MagForce1;
        u1Seq = [repmat(u1Sample1, 1, Td/2), repmat(u1Sample2, 1, Td/2)];
        u1IndPlot(SampleTimeStart:SampleTimeHalf) = u1IndSample1;
        u1IndPlot(SampleTimeHalf+1:SampleTimeEnd) = u1IndSample2;
    
        u2Sample1 = zeros(numCoil, 1); % Coil vector
        u2Sample2 = zeros(numCoil, 1);
        u2Sample1(u2IndSample1) = MagForce2;
        u2Sample2(u2IndSample2) = MagForce2;
    
        u2Seq = [repmat(u2Sample1, 1, Td/2), repmat(u2Sample2, 1, Td/2)];
        u2IndPlot(SampleTimeStart:SampleTimeHalf) = u2IndSample1;
        u2IndPlot(SampleTimeHalf+1:SampleTimeEnd) = u2IndSample2;
    
    
        uSeq = u1Seq + u2Seq;
    
    
        % Dynamics simulation loop
        for t = SampleTimeStart:SampleTimeEnd
            [x1Disk(t+1), y1Disk(t+1), vx1Disk(t+1), vy1Disk(t+1), F1Plot(1, t), F1Plot(2, t)] = ...
                UpdatePosVel(uSeq(:, t-SampleTimeStart+1), x1Disk(t), y1Disk(t), vx1Disk(t), vy1Disk(t), z1Disk, ...
                xGrid, yGrid, zGrid, numCoil, alpha, beta, dt);
            e1NormPlot(t) = norm(xy1RefSample - [x1Disk(t); y1Disk(t)]);
    
    
            [x2Disk(t+1), y2Disk(t+1), vx2Disk(t+1), vy2Disk(t+1), F2Plot(1, t), F2Plot(2, t)] = ...
                UpdatePosVel(uSeq(:, t-SampleTimeStart+1), x2Disk(t), y2Disk(t), vx2Disk(t), vy2Disk(t), z2Disk, ...
                xGrid, yGrid, zGrid, numCoil, alpha, beta, dt);
            e2NormPlot(t) = norm(xy2RefSample - [x2Disk(t); y2Disk(t)]);
    
            % end
        end

    else
        
        if (e1Sample<MagRng) 
            u1IndSample1 = u1RefIndSample;
            u1IndSample2 = u1RefIndSample;
        else
            [u1IndSample1, u1IndSample2] = SelectCoilFar(xy1Sample, u1RefIndSample, A1New, ...
                xyGridCell, numCoil, MagRng, s, e1Sample, xy1RefSample, 1);
        end
        if s==1
            [u2IndSample1, u2IndSample2] = SelectCoilFar(xy2Sample, u2RefIndSample, A2New, ...
                    xyGridCell, numCoil, MagRng, s, e2Sample, xy2Sample, 1);

            u2IndSample1_prev = u2IndSample1;
            u2IndSample2_prev = u2IndSample2;

            u2Sample1 = zeros(numCoil, 1); % Coil vector
        u2Sample2 = zeros(numCoil, 1);
        u2Sample1(u2IndSample1) = MagForce2;
        u2Sample2(u2IndSample2) = MagForce2;

        u2IndPlot(SampleTimeStart:SampleTimeHalf) = u2IndSample1;
        u2IndPlot(SampleTimeHalf+1:SampleTimeEnd) = u2IndSample2;
    
        u2Seq = [repmat(u2Sample1, 1, Td/2), repmat(u2Sample2, 1, Td/2)];
        

        else
            u2IndSample1 = u2IndSample2_prev;
            u2IndSample2 = u2IndSample2_prev;

            u2IndPlot(SampleTimeStart:SampleTimeHalf) = u2IndSample2;
            u2IndPlot(SampleTimeHalf+1:SampleTimeEnd) = u2IndSample2;
    
        u2Seq = [repmat(u2Sample2, 1, Td/2), repmat(u2Sample2, 1, Td/2)];
        end
        MagForce2 = 1;

        u1Sample1 = zeros(numCoil, 1); % Coil vector
        u1Sample2 = zeros(numCoil, 1);
        u1Sample1(u1IndSample1) = MagForce1;
        u1Sample2(u1IndSample2) = MagForce1;
        u1Seq = [repmat(u1Sample1, 1, Td/2), repmat(u1Sample2, 1, Td/2)];
        u1IndPlot(SampleTimeStart:SampleTimeHalf) = u1IndSample1;
        u1IndPlot(SampleTimeHalf+1:SampleTimeEnd) = u1IndSample2;

        uSeq = u1Seq + u2Seq;
    
    
        % Dynamics simulation loop
        for t = SampleTimeStart:SampleTimeEnd
            [x1Disk(t+1), y1Disk(t+1), vx1Disk(t+1), vy1Disk(t+1), F1Plot(1, t), F1Plot(2, t)] = ...
                UpdatePosVel(uSeq(:, t-SampleTimeStart+1), x1Disk(t), y1Disk(t), vx1Disk(t), vy1Disk(t), z1Disk, ...
                xGrid, yGrid, zGrid, numCoil, alpha, beta, dt);
            e1NormPlot(t) = norm(xy1RefSample - [x1Disk(t); y1Disk(t)]);
    
            
            [x2Disk(t+1), y2Disk(t+1), vx2Disk(t+1), vy2Disk(t+1), F2Plot(1, t), F2Plot(2, t)] = ...
            UpdatePosVel(uSeq(:, t-SampleTimeStart+1), x2Disk(t), y2Disk(t), vx2Disk(t), vy2Disk(t), z2Disk, ...
            xGrid, yGrid, zGrid, numCoil, alpha, beta, dt);
            e2NormPlot(t) = norm(xy2RefSample - [x2Disk(t); y2Disk(t)]);
    
            % end
        end

    end

end



% =========================================
% Plot Results
% =========================================

% Plot the trajectory
figure(99);
clf(99);
axis([[-1, 1]*(floor(GridSpacing*(xGridSize-1)/2)+1), ...
    [-1, 1]*(floor(GridSpacing*(yGridSize-1)/2)+1)]);
pbaspect([1, 1, 1]); daspect([1, 1, 1]);
box on; grid on;
hold on;
plot(xGrid(:), yGrid(:), 'g*');
plot(x1Disk, y1Disk, 'r-', 'LineWidth', 2);
plot(xy1RefPlot(1,:), xy1RefPlot(2,:), 'b--', 'LineWidth', 2);
plot(x2Disk, y2Disk, 'm-', 'LineWidth', 2);
plot(xy2RefPlot(1,:), xy2RefPlot(2,:), 'c--', 'LineWidth', 2);
title('Reference and Actual Trajectory of Moving Dipoles');
xlabel('X position');
ylabel('Y position');
legend('Coil grid', 'Disk 1 pos', 'Ref 1 pos', 'Disk 2 pos', 'Ref 2 pos', ...
    'location', 'northeastoutside');
hold off;

% Plot reference signals and errors
figure;
subplot(3, 2, 1);
hold on;
plot(vt, x1Disk, '-');
plot(vt, xy1RefPlot(1, :), '--');
title('X1 Reference and Actual Positions');
xlabel('Time Step');
ylabel('X Position');
hold off;

subplot(3, 2, 2);
hold on;
plot(vt, y1Disk, '-');
plot(vt, xy1RefPlot(2, :), '--');
title('Y1 Reference and Actual Positions');
xlabel('Time Step');
ylabel('Y Position');
hold off;

subplot(3, 2, 3);
plot(vt, vx1Disk, '-');
title('X1 Velocity');
xlabel('Time Step');
ylabel('X Velocity');

subplot(3, 2, 4);
plot(vt, vy1Disk, '-');
title('Y1 Velocity');
xlabel('Time Step');
ylabel('Y Velocity');

subplot(3, 2, 5);
hold on;
plot(vt, u1IndPlot, '-');
plot(vt, u1RefIndPlot, '--');
title('Reference and Actual Control Signals 1');
xlabel('Time Step');
ylabel('Control Signal');
hold off;

subplot(3, 2, 6);
plot(vt, e1NormPlot, '-');
title('Error Norm 1');
xlabel('Time Step');
ylabel('Error Norm');

figure;
subplot(3, 2, 1);
hold on;
plot(vt, x2Disk, '-');
plot(vt, xy2RefPlot(1, :), '--');
title('X2 Reference and Actual Positions');
xlabel('Time Step');
ylabel('X Position');
hold off;

subplot(3, 2, 2);
hold on;
plot(vt, y2Disk, '-');
plot(vt, xy2RefPlot(2, :), '--');
title('Y2 Reference and Actual Positions');
xlabel('Time Step');
ylabel('Y Position');
hold off;

subplot(3, 2, 3);
plot(vt, vx2Disk, '-');
title('X2 Velocity');
xlabel('Time Step');
ylabel('X Velocity');

subplot(3, 2, 4);
plot(vt, vy2Disk, '-');
title('Y2 Velocity');
xlabel('Time Step');
ylabel('Y Velocity');

subplot(3, 2, 5);
hold on;
plot(vt, u2IndPlot, '-');
plot(vt, u2RefIndPlot, '--');
title('Reference and Actual Control Signals 2');
xlabel('Time Step');
ylabel('Control Signal');
hold off;

subplot(3, 2, 6);
plot(vt, e2NormPlot, '-');
title('Error Norm 2');
xlabel('Time Step');
ylabel('Error Norm');

% Plot static dipoles and trajectory
figure;
axis([[-1, 1]*(floor(GridSpacing*(xGridSize-1)/2)+1), ...
    [-1, 1]*(floor(GridSpacing*(yGridSize-1)/2)+1)]);
pbaspect([1, 1, 1]); daspect([1, 1, 1]);
box on; grid on;
hold on;
plot(xGrid(:), yGrid(:), 'g*');

h1Disk = plot(x1Disk(1), y1Disk(1), 'ro', 'LineWidth', 2, 'MarkerSize', 12);  % Handle for actual plot
h1Ref = plot(xy1RefPlot(1,1), xy1RefPlot(2,1), 'b^', 'LineWidth', 2, 'MarkerSize', 12);  % Handle for reference plot
h1U = plot(xGrid(u1IndPlot(1)), yGrid(u1IndPlot(1)), 'rx', 'LineWidth', 2, 'MarkerSize', 12);  % Handle for active coil

h2Disk = plot(x2Disk(1), y2Disk(1), 'mo', 'LineWidth', 2, 'MarkerSize', 12);  % Handle for actual plot
h2Ref = plot(xy2RefPlot(1,1), xy2RefPlot(2,1), 'c^', 'LineWidth', 2, 'MarkerSize', 12);  % Handle for reference plot
h2U = plot(xGrid(u2IndPlot(1)), yGrid(u2IndPlot(1)), 'mx', 'LineWidth', 2, 'MarkerSize', 12);  % Handle for active coil

title('Reference and Actual Movement of Moving Dipoles');
xlabel('X position');
ylabel('Y position');
legend('Coil grid', 'Disk 1 pos', 'Ref 1 pos', 'Active coil 1', 'Disk 2 pos', 'Ref 2 pos', 'Active coil 2', ...
    'location', 'northeastoutside');

for k = 1:250:length(vt)-1
    set(h1Disk, 'XData', x1Disk(k), 'YData', y1Disk(k));  % Update actual plot
    set(h1Ref, 'XData', xy1RefPlot(1,k), 'YData', xy1RefPlot(2,k));  % Update reference plot
    set(h1U, 'XData', xGrid(u1IndPlot(k)), 'YData', yGrid(u1IndPlot(k))); % Update active coil

    set(h2Disk, 'XData', x2Disk(k), 'YData', y2Disk(k));  % Update actual plot
    set(h2Ref, 'XData', xy2RefPlot(1,k), 'YData', xy2RefPlot(2,k));  % Update reference plot
    set(h2U, 'XData', xGrid(u2IndPlot(k)), 'YData', yGrid(u2IndPlot(k))); % Update active coil

    drawnow;
    % Step-by-step simulation
%     if mod(k, 1000) == 1
%         disp(['Force = (', num2str(FPlot(1,k)), ',', num2str(FPlot(2,k)), ')']);
%         waitforbuttonpress;
%     end
end
hold off;

% =========================================
% Functions
% =========================================

function [uIndSample1, uIndSample2] = SelectCoilFar(xySample, uRefIndSample, ANew, ...
    xyGridCell, numCoil, MagRng, s, eSample, xyRefSample, diskInd, stop)
    % [Debug] Distance when out of range
    disp(['Disk ', num2str(diskInd), ' is out of range in loop ', num2str(s)]);
    disp(['Distance = ', num2str(eSample)]);
    if diskInd==1
        plot(xySample(1), xySample(2), 'ro', 'LineWidth', 2, 'MarkerSize', 12); 
        plot(xyRefSample(1), xyRefSample(2), 'b^', 'LineWidth', 2, 'MarkerSize', 12); 
    else
        plot(xySample(1), xySample(2), 'mo', 'LineWidth', 2, 'MarkerSize', 12); 
        plot(xyRefSample(1), xyRefSample(2), 'c^', 'LineWidth', 2, 'MarkerSize', 12); 
    end

    % Find the closest entry in the cell grid
    minDist = inf;
    minInd = 1;
    for i = 1:numCoil
        xyCell = xyGridCell{i};
        eCell = norm(xySample - xyCell);
        if eCell < minDist
            minDist = eCell;
            minInd = i;
        end
    end

    % Compute the shortest path using Dijkstra's algorithm
    
    [uSampleSeqCost, uSampleSeq] = dijkstra(ANew, minInd, uRefIndSample);
    % [Debug] Outputs from dijkstra
    disp(['Disk ', num2str(diskInd), ' from ', num2str(minInd), ' to ', num2str(uRefIndSample)]);
    disp(['dijkstra seq = ', num2str(uSampleSeq)]);
    disp(['dijkstra cost = ', num2str(uSampleSeqCost)]);
    uSampleSeq = fliplr(uSampleSeq);
    % To-do: choose the best in the closest four

    % Choose the 2nd and 3rd coils if close enough
    if norm(xySample - xyGridCell{uSampleSeq(2)}) <= MagRng
        uIndSample1 = uSampleSeq(2);
        uIndSample2 = uSampleSeq(3);
    else % Else choose the 1st and 2nd coils
        uIndSample1 = uSampleSeq(1);
        uIndSample2 = uSampleSeq(2);
    end
end

function ANew = AdjustWeight(A, xyPos, weight, xyGridCell, xGridSize, yGridSize, MagRng)
% Change the weight of edges to coils within MagRng of xyPos
    ANew = A;
    for i = 1:xGridSize
        for j = 1:yGridSize
            xyCell = xyGridCell{i, j};
            if norm(xyPos - xyCell) <= MagRng
                currentInd = sub2ind([xGridSize, yGridSize], i, j);
                neighbors = [i, j-1; i, j+1; i-1, j; i+1, j; ...
                    i-1, j-1; i-1, j+1; i+1, j-1; i+1, j+1];
                for k = 1:length(neighbors)
                    ni = neighbors(k, 1);
                    nj = neighbors(k, 2);
                    if ni >= 1 && ni <= xGridSize && nj >= 1 && nj <= yGridSize
                        neighborInd = sub2ind([xGridSize, yGridSize], ni, nj);
                        ANew(currentInd, neighborInd) = weight;
                        ANew(neighborInd, currentInd) = weight;
                    end
                end
            end
        end
    end
end

function [xDiskNew, yDiskNew, vxDiskNew, vyDiskNew, Fx_total, Fy_total] = ...
    UpdatePosVel(uSeq, xDisk, yDisk, vxDisk, vyDisk, zDisk, ...
    xGrid, yGrid, zGrid, numCoil, alpha, beta, dt)
    % Calculate forces
    Fx_total = 0;
    Fy_total = 0;
    for i = 1:numCoil
        % [To-do] Should MagRng be enfored in dynamics simulation?
        % if norm([xDisk, yDisk] - [xGrid(i), yGrid(i)]) <= MagRng
            r = norm([xDisk, yDisk, zDisk] - [xGrid(i), yGrid(i), zGrid(i)]); % Euclidean distance
            F = alpha / r^5 * (1 - 5 * (zDisk - zGrid(i))^2 / r^2);
    
            Fx_total = Fx_total + uSeq(i) * F * (xDisk - xGrid(i));  % Sum of forces in x direction
            Fy_total = Fy_total + uSeq(i) * F * (yDisk - yGrid(i));  % Sum of forces in y direction
        % end
    end

    % Update position and velocity
    xDiskNew = xDisk + vxDisk * dt;
    yDiskNew = yDisk + vyDisk * dt;
    vxDiskNew = vxDisk - beta * vxDisk * dt + Fx_total * dt;  % Force component in x
    vyDiskNew = vyDisk - beta * vyDisk * dt + Fy_total * dt;  % Force component in y
end
