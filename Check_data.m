%load('lorenz_nonst.mat');

%Check the lengths below
for i=1:size(results,1)
    a = results{i};
    
    [unique_vals, ~, idx] = unique(a(:,7)); 
    counts = histc(a(:,7), unique_vals)
end



% For 1:0.05:500
    % 9,981 per sim
    % 9,981 per sim per rho * 100 simulations
    % multiply by number of rhos

% For 1:0.05:5000
    % 9998100 points, 100 simulations
    % 99,981 points per sim per rho  * simulation
    % Multiple by number of rhos

% Detailed:
    % inc/dec - 11 rhos (at the top not doubled)
        % Rhos = 25 65  105 145 185 225 185 145 105 65 25
    % not inc - 6 rhos (increasing from 25 225 last)
        % Rhos = 25    65   105   145   185   225 (6)
    %~9090 points per rho value in inc/dec
    %~16665 points per rho value in inc

% inc/dec - 11 rhos (at the top not doubled)
% not inc - 6 rhos (increasing from 25 225 last)

%~ 1817841 / 2,  908922 / 2 points per rho value in inc/dec
%~ 1666354 points per rho value in inc
% In total, divide by 100 for each sim.
% 18178, 9089 for inc dec
% 16663 for not inc