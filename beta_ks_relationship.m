clc; clear;
% This program finds the relationship between beta, n, and k in Power-law fluids
% as the result from eqn 1.327, p.127 Non-Newtonian Flow and Applied Rheology Engineering by
% R. P. Chhabra, J. F. Richardson -- Butterworth-Heinemann_IChemE, 2ed, 2008

%% Parameters

% Read k, n, and beta values from the CSV file
data = readtable('results.csv');

% Initialize figure
figure;
hold on;

% Define markers and colors for different k values
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*', 'x'};
colors = lines(length(unique(data.k))); % Use distinct colors for each unique k value

% Initialize a vector to store all ks values
all_ks_values = [];

unique_k_values = unique(data.k);
for i = 1:length(unique_k_values)
    k = unique_k_values(i);
    k_filter = data.k == k;
    
    n_values = data.n(k_filter);
    beta_values = data.beta(k_filter);
    ks_values = [];

    for j = 1:length(beta_values)
        beta = beta_values(j);
        n = n_values(j);

        beta1_sq = (k^2-1)/(2*log(k)); % Given k
        beta1 = sqrt(beta1_sq);     % Fixed parameter beta1
        s = 1 / n; % Define s = 1/n

        % Compute the first term
        term1 = (1 + k^2 - 2 * beta1^2)/4;
        part1 = (term1)^(1/(n-1));
        % Compute the second term
        term2_ts = (3+s)*(1-k^2);
        term2_ms_p1 = (1-beta^2)^(1+s);
        term2_ms_p2 = k^(1-s);
        term2_ms_p3 = (beta^2-k^2)^(1+s);

        term2_ms = term2_ms_p1 - term2_ms_p2*term2_ms_p3;

        term2 = term2_ts/term2_ms;

        part2 = term2^(n/(n-1));

        % Compute Ks
        ks = (1-k) * part1 * part2;
        if ks >= 0 % Only keep non-negative Ks values
            ks_values = [ks_values, ks];
            all_ks_values = [all_ks_values, ks]; % Store ks in the all_ks_values vector
        else
            beta_values(j) = NaN; % Mark beta as NaN to remove it later
        end
    end

    % Remove NaN values from beta_values
    beta_values = beta_values(~isnan(beta_values));
    
    % Plot ks values for the current k as a line with filled markers
    plot( ks_values,beta_values, '-o', ...
        'Marker', markers{mod(i-1, length(markers))+1}, ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', colors(i,:), ...
        'Color', colors(i,:), ...
        'DisplayName', sprintf('k = %.2f', k));
end

% Calculate the mean of all ks values
mean_ks = mean(all_ks_values)
disp(['Mean Ks: ', num2str(mean_ks)]);

% Set plot properties
ylabel('\beta');
xlabel('k_s');
title('Relationship between \beta and k_s for different values of k and n');
legend show;
grid on;
hold off;