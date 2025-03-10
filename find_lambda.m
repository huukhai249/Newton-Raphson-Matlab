clc; clear;
% This program finds the relationship between lambda, n, and sigma in Power-law fluids
% as the result from eqn 1.327, p.127 Non-Newtonian Flow and Applied Rheology Engineering by 
% R. P. Chhabra, J. F. Richardson -- Butterworth-Heinemann_IChemE, 2ed, 2008 

%% Parameters

lambda_start = 0.1;
lambda_end = 1;
sigma_start = 0.1;
sigma_end = 0.9;
n_start =0.1;
n_end = 2.0;

sigma_range = sigma_start:0.05:sigma_end; % Range for sigma
n_range = n_start:0.1:n_end; % Range for n
lambda_range = lambda_start:0.1:lambda_end; % Range for lambda (initial guess)

% Tolerances and iteration limit for Newton-Raphson
tolerance = 1e-16;
max_iterations = 100;
h = 1e-6; % Step size for numerical derivative

% Initialize arrays to store results
n_values = [];
lambda_values = [];
sigma_values = [];

% Loop over all combinations of sigma, lambda, and n
for sigma = sigma_range;
    for n = n_range
        best_lambda = NaN;
        best_error = inf;
        for lambda_init = lambda_range
            lambda = lambda_init; % Initial guess for lambda
            
            % Newton-Raphson iteration
            for iter = 1:max_iterations
                % Define the function f(lambda) = I1 - I2
                % Numerical integration for I1: int_sigma^lambda (lambda^2/x^2 - x)^(1/n) dx
                x1 = linspace(sigma, lambda, 200); % Discretize the interval
                integrand1 = ((lambda^2 ./ x1)- x1).^(1/n);
                I1 = trapz(x1, integrand1);
                
                % Numerical integration for I2: int_lambda^1 (x - lambda^2/x)^(1/n) dx
                x2 = linspace(lambda, 1, 200); % Discretize the interval
                integrand2 = (x2 - (lambda^2 ./ x2)).^(1/n);
                I2 = trapz(x2, integrand2);
                
                % Function f(lambda)
                f = I1 - I2;
                
                % Numerical derivative f'(lambda) using central difference
                lambda_plus_h = lambda + h;
                lambda_minus_h = lambda - h;
                
                % Compute f(lambda + h)
                % f(lambda + h) = I1_plus -I2_plus

                x1_plus = linspace(sigma, lambda_plus_h, 100);
                integrand1_plus = (lambda_plus_h^2 ./ x1_plus - x1_plus).^(1/n);
                I1_plus = trapz(x1_plus, integrand1_plus);
                
                x2_plus = linspace(lambda_plus_h, 1, 100);
                integrand2_plus = (x2_plus - lambda_plus_h^2 ./ x2_plus).^(1/n);
                I2_plus = trapz(x2_plus, integrand2_plus);
                f_plus = I1_plus - I2_plus;
                
                % Compute f(lambda - h)
                x1_minus = linspace(sigma, lambda_minus_h, 100);
                integrand1_minus = (lambda_minus_h^2 ./ x1_minus - x1_minus).^(1/n);
                I1_minus = trapz(x1_minus, integrand1_minus);
                
                x2_minus = linspace(lambda_minus_h, 1, 100);
                integrand2_minus = (x2_minus - lambda_minus_h^2 ./ x2_minus).^(1/n);
                I2_minus = trapz(x2_minus, integrand2_minus);
                f_minus = I1_minus - I2_minus;
                
                % Central difference approximation
                df = (f_plus - f_minus) / (2 * h);
                
                % Check for division by zero or NaN
                if abs(df) < 1e-10
                    disp('Derivative too small, breaking loop.');
                    break;
                end
                
                %%%%%%%%%%%%%%%%%%% Update lambda
                lambda_new = lambda - f / df;
                %%%%%%%%% This  is Newton-Rapshon interator %%%%%%%%%%%%
                
                % Check convergence
                if abs(lambda_new - lambda) < tolerance
                    % Check if this lambda is the best one for this n
                    error = abs(f);
                    if error < best_error
                        best_error = error;
                        best_lambda = lambda_new;
                    end
                    break;
                end
                
                lambda = lambda_new;
                
                % Check if lambda is within bounds (0.1 to 0.9)
                if lambda < lambda_start || lambda > lambda_end
                    disp('Lambda out of bounds, restarting with new initial guess.');
                    break;
                end
                
                if iter == max_iterations
                    disp('Max iterations reached, no convergence.');
                end
            end
        end
        
        % Store the best lambda for this n
        if ~isnan(best_lambda)
            n_values = [n_values, n];
            lambda_values = [lambda_values, best_lambda];
            sigma_values = [sigma_values, sigma];
            fprintf('Best result for sigma = %.2f, n = %.2f: lambda = %.6f\n', sigma, n, best_lambda);
        end
    end
end

% Create the plot
figure; % Create a new figure
hold on; % Allow multiple plots on the same figure

% Define markers or colors for different sigma values
markers = {'o', 's', '^', 'd', 'p', 'h', '*', 'x'}; % Different markers for each sigma
colors = lines(length(sigma_range)); % Use distinct colors

% Plot for each unique sigma value
for i = 1:length(sigma_range)
    sigma_idx = sigma_values == sigma_range(i);
    plot(n_values(sigma_idx), lambda_values(sigma_idx), ...
         'Marker', markers{mod(i-1, length(markers))+1}, ...
         'LineStyle', 'none', ...
         'Color', colors(i,:), ...
         'DisplayName', sprintf('k = %.2f', sigma_range(i)));
end

% Add labels and legend
xlabel('n');
ylabel('\beta');
title('Relationship between n and \beta for different k values');
legend('show');
grid on;

% Ensure the plot is visually appealing
axis([n_start n_end lambda_start lambda_end]); % Set axis limits based on the ranges (0.1 to 0.9)
hold off;