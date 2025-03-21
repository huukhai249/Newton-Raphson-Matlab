clc; clear;
% This program finds the relationship between beta, n, and k in Power-law fluids
% as the result from eqn 1.327, p.127 Non-Newtonian Flow and Applied Rheology Engineering by 
% R. P. Chhabra, J. F. Richardson -- Butterworth-Heinemann_IChemE, 2ed, 2008 

%% Parameters

beta_start = 0.1;
beta_end = 1;
k_start = 0.1;
k_end = 0.9;
n_start = 0.1;
n_end = 1.0;

k_range = k_start:0.05:k_end; % Range for k
n_range = n_start:0.05:n_end; % Range for n
beta_range = beta_start:0.1:beta_end; % Range for beta (initial guess)

% Tolerances and iteration limit for Newton-Raphson
tolerance = 1e-12;
max_iterations = 100;
h = 1e-6; % Step size for numerical derivative

% Initialize arrays to store results
n_values = [];
beta_values = [];
k_values = [];

% Loop over all combinations of k, beta, and n
for k = k_range
    for n = n_range
        best_beta = NaN;
        best_error = inf;
        for beta_init = beta_range
            beta = beta_init; % Initial guess for beta
            
            % Newton-Raphson iteration
            for iter = 1:max_iterations
                % Define the function f(beta) = I1 - I2
                % Numerical integration for I1: int_k^beta (beta^2/x^2 - x)^(1/n) dx
                x1 = linspace(k, beta, 200); % Discretize the interval
                integrand1 = ((beta^2 ./ x1) - x1).^(1/n);
                I1 = trapz(x1, integrand1);
                
                % Numerical integration for I2: int_beta^1 (x - beta^2/x)^(1/n) dx
                x2 = linspace(beta, 1, 200); % Discretize the interval
                integrand2 = (x2 - (beta^2 ./ x2)).^(1/n);
                I2 = trapz(x2, integrand2);
                
                % Function f(beta)
                f = I1 - I2;
                
                % Numerical derivative f'(beta) using central difference
                beta_plus_h = beta + h;
                beta_minus_h = beta - h;
                
                % Compute f(beta + h)
                % f(beta + h) = I1_plus - I2_plus

                x1_plus = linspace(k, beta_plus_h, 100);
                integrand1_plus = (beta_plus_h^2 ./ x1_plus - x1_plus).^(1/n);
                I1_plus = trapz(x1_plus, integrand1_plus);
                
                x2_plus = linspace(beta_plus_h, 1, 100);
                integrand2_plus = (x2_plus - beta_plus_h^2 ./ x2_plus).^(1/n);
                I2_plus = trapz(x2_plus, integrand2_plus);
                f_plus = I1_plus - I2_plus;
                
                % Compute f(beta - h)
                x1_minus = linspace(k, beta_minus_h, 100);
                integrand1_minus = (beta_minus_h^2 ./ x1_minus - x1_minus).^(1/n);
                I1_minus = trapz(x1_minus, integrand1_minus);
                
                x2_minus = linspace(beta_minus_h, 1, 100);
                integrand2_minus = (x2_minus - beta_minus_h^2 ./ x2_minus).^(1/n);
                I2_minus = trapz(x2_minus, integrand2_minus);
                f_minus = I1_minus - I2_minus;
                
                % Central difference approximation
                df = (f_plus - f_minus) / (2 * h);
                
                % Check for division by zero or NaN
                if abs(df) < 1e-10
                    disp('Derivative too small, breaking loop.');
                    break;
                end
                
                % Update beta
                beta_new = beta - f / df;
                % This is the Newton-Rapshon iterator
                
                % Check convergence
                if abs(beta_new - beta) < tolerance
                    % Check if this beta is the best one for this n
                    error = abs(f);
                    if error < best_error
                        best_error = error;
                        best_beta = beta_new;
                    end
                    break;
                end
                
                beta = beta_new;
                
                % Check if beta is within bounds (0.1 to 0.9)
                if beta < beta_start || beta > beta_end
                    disp('Beta out of bounds, restarting with new initial guess.');
                    break;
                end
                
                if iter == max_iterations
                    disp('Max iterations reached, no convergence.');
                end
            end
        end
        
        % Store the best beta for this n
        if ~isnan(best_beta)
            n_values = [n_values, n];
            beta_values = [beta_values, best_beta];
            k_values = [k_values, k];
            fprintf('Best result for k = %.2f, n = %.2f: beta = %.6f\n', k, n, best_beta);
        end
    end
end

% Create the plot for n and beta
figure; % Create a new figure
hold on; % Allow multiple plots on the same figure

% Define markers or colors for different k values
markers = {'o', 's', '^', 'd', 'p', 'h', '*', 'x'}; % Different markers for each k
colors = lines(length(k_range)); % Use distinct colors

% Plot for each unique k value
for i = 1:length(k_range)
    k_idx = k_values == k_range(i);
    plot(n_values(k_idx), beta_values(k_idx), ...
         'Marker', markers{mod(i-1, length(markers))+1}, ...
         'LineStyle', 'none', ...
         'Color', colors(i,:), ...
         'DisplayName', sprintf('k = %.2f', k_range(i)));
end

% Add labels and legend
xlabel('n');
ylabel('\beta');
title('Relationship between n and \beta for different k values');
legend('show');
grid on;

% Ensure the plot is visually appealing
axis([n_start n_end beta_start beta_end]); % Set axis limits based on the ranges (0.1 to 0.9)
hold off;

% Create the scatter plot for unique_k_values and mean_beta_values
figure; % Create a new figure
scatter3(k_values, n_values, beta_values, 'filled');
xlabel('k');
ylabel('n');
zlabel('\beta');
title('Scatter Plot of \beta vs. k and n');
grid on;

% Perform surface fitting
fit_type = fittype('poly11'); % Linear polynomial in both k and n
fit_result = fit([k_values(:), n_values(:)], beta_values(:), fit_type);

% Generate fitting values
[k_fit, n_fit] = meshgrid(linspace(min(k_values), max(k_values), 100), linspace(min(n_values), max(n_values), 100));
beta_fit = fit_result.p00 + fit_result.p10 * k_fit + fit_result.p01 * n_fit;

% Plot the fitting surface
hold on;
mesh(k_fit, n_fit, beta_fit, 'EdgeColor', 'red', 'DisplayName', 'Fitting Surface');
legend('show');
hold off;

% Display the fitting result
disp('Fitting result:');
disp(fit_result);

% Display the fitting polynomial expression
fprintf('Fitting Polynomial: beta(k, n) = %.4f + %.4f*k + %.4f*n\n', fit_result.p00, fit_result.p10, fit_result.p01);