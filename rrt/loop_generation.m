clear all
close all
clc

loop = generate_valid_loop(15, 50);  % 15x15 grid, 50-step loop



function loop_indices = generate_valid_loop(gridSize, loopLength, startIdx)
    if nargin < 3
        startIdx = randi(gridSize^2);
    end
    
    visited = false(gridSize);
    [r, c] = ind2sub([gridSize, gridSize], startIdx);
    visited(r, c) = true;
    loop = startIdx;

    for i = 2:loopLength
        neighbors = get_neighbors_8(r, c, gridSize);
        valid = neighbors(~visited(sub2ind([gridSize, gridSize], neighbors(:,1), neighbors(:,2))), :);
        if isempty(valid)
            break;
        end
        next = valid(randi(size(valid,1)), :);
        visited(next(1), next(2)) = true;
        loop(end+1) = sub2ind([gridSize, gridSize], next(1), next(2));
        r = next(1); c = next(2);
    end

    % Close the loop if possible
    [sr, sc] = ind2sub([gridSize, gridSize], loop(1));
    if max(abs(r - sr), abs(c - sc)) == 1
        loop(end+1) = loop(1);
    end

    loop_indices = loop;
end

function neighbors = get_neighbors_8(r, c, gridSize)
    directions = [-1 -1; -1 0; -1 1;
                   0 -1;         0 1;
                   1 -1;  1 0;  1 1];
    neighbors = [];
    for i = 1:size(directions,1)
        nr = r + directions(i,1);
        nc = c + directions(i,2);
        if nr >= 1 && nr <= gridSize && nc >= 1 && nc <= gridSize
            neighbors(end+1, :) = [nr, nc];
        end
    end
end
