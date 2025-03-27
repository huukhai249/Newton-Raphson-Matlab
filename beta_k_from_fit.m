clc; clear;
% This script calculates and plots ks values for beta for k values from 0.1 to 0.9
% and n values from 0.1 to 0.9.

% Parameters
k_start = 0.1;
k_end = 0.9;
k_step = 0.1;
k_values = k_start:k_step:k_end;
n_values = 0.1:0.1:0.9; % n values from 0.1 to 0.9

% Coefficients from the given polynomial fit for beta
beta_values = zeros(size(k_values));

% Initialize figure
figure;
hold on;

% Define markers and colors for different n values
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*', 'x'};
colors = lines(length(n_values)); % Use distinct colors for each unique n value

% Loop over each n value
for j = 1:length(k_values)
    ks_values = zeros(size(n_values));
    k = k_values(j);

    % Loop over each k value to calculate ks
    for i = 1:length(n_values)
        n = n_values(i);
        beta = 0.3410*k^2 + 0.9106*k - 0.2683;
        beta_values(i)=beta;
        % Skip negative beta values
        if beta <= 0
            continue;
        end

        beta1_sq = (k^2 - 1) / (2 * log(k)); % Given k
        beta1 = sqrt(beta1_sq); % Fixed parameter beta1
        s = 1 / n; % Define s = 1/n

        % Compute the first term
        term1 = (1 + k^2 - 2 * beta1^2) / 4;
        part1 = (term1)^(1 / (n - 1));

        % Compute the second term
        term2_ts = (3 + s) * (1 - k^2);
        term2_ms_p1 = (1 - beta^2)^(1 + s);
        term2_ms_p2 = k^(1 - s);
        term2_ms_p3 = (beta^2 - k^2)^(1 + s);
        term2_ms = term2_ms_p1 - term2_ms_p2 * term2_ms_p3;
        term2 = term2_ts / term2_ms;

        part2 = term2^(n / (n - 1));

        % Compute ks
        ks = (1 - k) * part1 * part2;

        % Skip negative ks values
        if ks <= 0
            continue;
        end

        ks_values(i) = ks;
    end
    
    % Remove zeros from ks_values and corresponding beta_values
    valid_idx = ks_values > 0;
    ks_values = ks_values(valid_idx);
    beta_values_valid = beta_values(valid_idx);

    % Plot ks values for the current n as a line with filled markers
    plot( ks_values, beta_values_valid, '-o', ...
        'Marker', markers{mod(j-1, length(markers))+1}, ...
        'MarkerFaceColor', colors(j,:), ...
        'MarkerEdgeColor', colors(j,:), ...
        'Color', colors(j,:), ...
        'DisplayName', sprintf('k = %.1f', k));
end

% Set plot properties
ylabel('\beta');
xlabel('k_s');
title('k_s vs. \beta for different values of n');
legend show;
set(gca, 'XScale', 'log');

grid on;
hold off;