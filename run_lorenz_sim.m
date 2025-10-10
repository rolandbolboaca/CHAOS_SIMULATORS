clear; clc;

% Default Params
sigma = 10;
beta = 8/3;
rho = 25;

rho_range = [25 25 225];

% Time span with delta_t
tspan = 1:0.1:1000;

% Initial conditions
xinit = [1 1 1];

rho_non_st = 0; %1 or 2 or 0
step_rho = 0;

% 1 Normal, 2 Uniform
noise = 1;

% IC Ranges
RANGE = [-5 5];

x_range = RANGE;
y_range = RANGE;
z_range = RANGE;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generate data with stepwise rho values
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotHOSandHists(s, res);

sims = 100;
tspan = 1:0.005:100;

rho_range = [25 40 250];
% rho_range = [150 150 150];
% % Params: Increase-Decrease, Random shift rho
% % ID: 20,30,40,30,20, etc.
% % RS: 40, 20, 30, etc.
% % ID and RS: more simulations and more rhos.
% % false, true will generate random shifts
% paramList = {
%     true, true;
%     true,  false;
%     false, true;
%     false, false
% };

RANDOM_LENGTH = true;
paramList = {
     false, true;
};

numRuns = size(paramList,1) * 2;   % 2 runs for each config
results = cell(numRuns,1);

parfor k = 1:numRuns

    inc_dec = paramList{ceil(k/2),1};
    random_shuffle = paramList{ceil(k/2),2};

    [s, res] = generateSimulations(sigma, beta, rho, tspan, xinit, noise, step_rho, ...
                                   sims, x_range, y_range, x_range, rho_range, ...
                                   inc_dec, random_shuffle, RANDOM_LENGTH);
    results{k} = res;
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

% [(max - min)/step] + 1
% 0.005 = 19801 with 6 rho = ~ 
save('lorenz_nonst.mat', '-v7.3'); 

function [s, res] = generateSimulations(sigma, beta, rho, tspan, xinit, noise, step_rho, ...
                                        sims, x_range, y_range, z_range, rho_range, inc_dec, ...
                                        random_shuffle, RANDOM_LENGTH)

    s = LorenzSystemSim(sigma, beta, rho, tspan, xinit, noise, step_rho);
    res = s.run_extended_random_ic_rho_intervals_cont(sims, x_range, y_range, z_range, ... 
                                                        rho_range, inc_dec, random_shuffle, ...
                                                        RANDOM_LENGTH);
    
    name = createName(noise, rho_range, inc_dec, x_range, sims, random_shuffle);
    s.writeToCSVFile(res, name);

end

function out = createName(noise, rho_range, inc_dec, RANGE, sims, random_shuffle)

    if noise == 1
        noise_n = "G";
    else
        noise_n = "U";
    end

    rho_n = "Rho_";
    rho_n = rho_n + string(rho_range(1)) + "_" + string(rho_range(2)) + "_" + string(rho_range(3));
    
    inc_dec_n = "";
    if inc_dec == true
        inc_dec_n = "ID";
    end
    date_n = string(datetime("now", "Format",'dd_MM_yyyy_HH_mm_SS'));

    rs = "";
    if random_shuffle == true
        rs = "RS";
    end

    out = "Norm_IC_" + string(RANGE(1)) + string(RANGE(2)) + ...
          "_" + noise_n + "_" + rho_n + "_" + inc_dec_n + ...
          "_" + rs + "_" + "Sim_" + string(sims) + "_" + ... 
          date_n + ".csv";
end

function plotHOSandHists(s, res)

    seg = s.getSegmenLength;

    % % Calculate the number of segments
    num_segments = round(size(res, 1) / seg);

    % Define the column indices for means, std, kurtosis, and skewness
    columns = 4:6;

    % Preallocate arrays for means, stds, kurtosis, and skewness
    means = zeros(num_segments, length(columns));
    stands = zeros(num_segments, length(columns));
    kurts = zeros(num_segments, length(columns));
    skews = zeros(num_segments, length(columns));

    idx = 1;
    for i = 1:seg:size(res,1) - seg

        start = i;
        end_r = i + seg;
        current_range = start:end_r; 

        means(idx, :) = mean(res(current_range, columns));
        stands(idx, :) = std(res(current_range, columns));
        skews(idx, :) = skewness(res(current_range, columns));
        kurts(idx, :) = kurtosis(res(current_range, columns));

        idx = idx + 1;
    end

    % 
    % figure; plot(means(1:end-1,:)); legend(["x","y","z"]);title("Mean Values");
    % figure; plot(stands(1:end-1,:)); legend(["x","y","z"]);title("Standard Deviataion Values");
    % figure; plot(skews(1:end-1,:)); legend(["x","y","z"]);title("Skewness Values");
    % figure; plot(kurts(1:end-1,:)); legend(["x","y","z"]);title("Kurtosis Values");
    figure; plot(res(:,4:end));

    % st = 8;
    % en = 9;
    % 
    % figure;histogram(res(seg*st:seg*en,4), 100)
    % figure;histogram(res(seg*st:seg*en,5), 100)
    % figure;histogram(res(seg*st:seg*en,6), 100)

end