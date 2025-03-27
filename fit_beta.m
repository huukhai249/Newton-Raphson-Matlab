clc; clear;
% This script reads k, n, and beta values from results.csv, calculates the min, max, and mean
% of beta for each unique k, and fits beta as a function of k.

% Read k, n, and beta values from the CSV file
data = readtable('results.csv');

% Get unique k values
unique_k_values = unique(data.k);

% Initialize arrays to store min, max, and mean beta values
min_beta_values = zeros(size(unique_k_values));
max_beta_values = zeros(size(unique_k_values));
mean_beta_values = zeros(size(unique_k_values));

% Loop over each unique k value to calculate min, max, and mean beta
for i = 1:length(unique_k_values)
    k = unique_k_values(i);
    k_filter = data.k == k;
    
    % Extract beta values for the current k
    beta_values = data.beta(k_filter);
    
    % Calculate the min, max, and mean beta for the current k
    min_beta_values(i) = min(beta_values);
    max_beta_values(i) = max(beta_values);
    mean_beta_values(i) = mean(beta_values);
end

% Perform fitting: beta = f(k)
[fit_result, gof] = fit(unique_k_values, mean_beta_values, 'poly2'); % Quadratic fit

% Display the fitting result
disp('Fitting result:');
disp(fit_result);

% Display the goodness of fit
disp('Goodness of fit:');
disp(gof);

% Display the fitting polynomial expression
fprintf('Fitting Polynomial: beta(k) = %.4f + %.4f*k + %.4f*k^2\n', fit_result.p1, fit_result.p2, fit_result.p3);

% Plot the min, max, and mean beta values and the fitting curve
figure;
hold on;

% Plot the min and max beta values as shaded area
fill([unique_k_values; flipud(unique_k_values)], [min_beta_values; flipud(max_beta_values)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'DisplayName', 'Min-Max Range');

% Plot the mean beta values
scatter(unique_k_values, mean_beta_values, 'filled', 'DisplayName', 'Mean \beta values');

% Plot the fitting curve
plot(unique_k_values, fit_result(unique_k_values), 'r-', 'LineWidth', 2, 'DisplayName', 'Fitting curve');

xlabel('k');
ylabel('\beta');
title('Min, Max, Mean \beta vs. k and Fitting Curve');
legend show;
grid on;
hold off;