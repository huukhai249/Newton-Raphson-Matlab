clc; clear;

%% Parameters
lambda_start = 0.1; % in Xuesi's work, lambda is beta
lambda_end = 1;
n_start = 0.1;
n_end = 2.0;
sigma_fixed = 0.36881; % Fixed sigma value (or k value in Xuesi's work)
% sigma_fixed take from find_sigma.m program
%% ----------  START NEWTON-RAPHSON LOOP  --------------%%

n_range = n_start:0.01:n_end; % Range for n
lambda_range = lambda_start:0.01:lambda_end; % Range for lambda (initial guess)

% Tolerances and iteration limit for Newton-Raphson
tolerance = 1e-12; % Adjusted tolerance
max_iterations = 100;
h = 1e-6; % Step size for numerical derivative

% Initialize arrays to store results
n_values = [];
lambda_values = [];
f_values = [];

fprintf('Iteration\tSigma\t\tN\t\tLambda Init\tLambda\n');
fprintf('-------------------------------------------------------------\n');

%% Loop over all combinations of lambda and n
for n = n_range
    found = false;
    for lambda_init = lambda_range
        lambda = lambda_init; % Initial guess for lambda
        
        % Newton-Raphson iteration
        for iter = 1:max_iterations
            % Define the function f(lambda) = I1 - I2
            % Numerical integration for I1: int_sigma^lambda (lambda^2/x^2 - x)^(1/n) dx
            x1 = linspace(sigma_fixed, lambda, 200); % Discretize the interval
            integrand1 = ((lambda^2 ./ x1) - x1).^(1/n);
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
            x1_plus = linspace(sigma_fixed, lambda_plus_h, 100);
            integrand1_plus = (lambda_plus_h^2 ./ x1_plus - x1_plus).^(1/n);
            I1_plus = trapz(x1_plus, integrand1_plus);
            
            x2_plus = linspace(lambda_plus_h, 1, 100);
            integrand2_plus = (x2_plus - lambda_plus_h^2 ./ x2_plus).^(1/n);
            I2_plus = trapz(x2_plus, integrand2_plus);
            f_plus = I1_plus - I2_plus;
            
            % Compute f(lambda - h)
            x1_minus = linspace(sigma_fixed, lambda_minus_h, 100);
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
            
            % Update lambda
            lambda_new = lambda - f / df;
            
            % Display current guess and function value
            fprintf('%d\t\t%.6f\t%.2f\t%.2f\t%.6f\n', iter, sigma_fixed, n, lambda_init, lambda_new);
            
            % Check convergence
            if abs(lambda_new - lambda) < tolerance
                % Store the converged values
                n_values = [n_values, n];
                lambda_values = [lambda_values, lambda_new];
                f_values = [f_values, abs(f)]; % Store the absolute value of f(lambda)
                found = true;
                break;
            end
            
            lambda = lambda_new;
            
            % Check if lambda is within bounds (0.1 to 1)
            if lambda < lambda_start || lambda > lambda_end
                disp('Lambda out of bounds, restarting with new initial guess.');
                break;
            end
            
            if iter == max_iterations
                disp('Max iterations reached, no convergence.');
            end
        end
        if found
            break;
        end
    end
end

%% Filter the results to get one lambda value for each n
filtered_n_values = [];
filtered_lambda_values = [];

unique_n_values = unique(n_values);
for i = 1:length(unique_n_values)
    n = unique_n_values(i);
    idx = find(n_values == n);
    if length(idx) > 1
        [~, min_idx] = min(f_values(idx)); % Find the index with the minimum f value
        best_idx = idx(min_idx);
    else
        best_idx = idx;
    end
    filtered_n_values = [filtered_n_values, n_values(best_idx)];
    filtered_lambda_values = [filtered_lambda_values, lambda_values(best_idx)];
end

%% Separate the results into two groups: n < 1 and n >= 1
n_less_than_1 = filtered_n_values(filtered_n_values < 1);
lambda_less_than_1 = filtered_lambda_values(filtered_n_values < 1);

n_greater_than_1 = filtered_n_values(filtered_n_values >= 1);
lambda_greater_than_1 = filtered_lambda_values(filtered_n_values >= 1);

%% Create the plot for n < 1
figure; % Create a new figure
hold on; % Allow multiple plots on the same figure

% Plot the results for n < 1
plot(n_less_than_1, lambda_less_than_1, 'bo');
xlabel('n');
ylabel('\beta');
title('Relationship between n and \beta for n < 1, fixed k=0.368694, Newtonian fluid');
grid on;

% Ensure the plot is visually appealing
axis([n_start 1 lambda_start lambda_end]); % Set axis limits based on the ranges (0.1 to 1)
hold off;

% Write results to an Excel file for n < 1
filename_less_than_1 = 'lambda_n_lt1.xlsx';

if isfile(filename_less_than_1)
    delete(filename_less_than_1); % Delete the existing file if it exists
end

result_table_less_than_1 = table(lambda_less_than_1', n_less_than_1', 'VariableNames', {'Lambda', 'N'});
writetable(result_table_less_than_1, filename_less_than_1);

%% Create the plot for n >= 1
% figure; % Create a new figure
% hold on; % Allow multiple plots on the same figure
% 
% % Plot the results for n >= 1
% plot(n_greater_than_1, lambda_greater_than_1, 'ro-');
% xlabel('n');
% ylabel('\beta');
% title('Relationship between n and \beta for n >= 1, fixed k=0.368694, Newtonian fluid');
% grid on;
% 
% % Ensure the plot is visually appealing
% axis([1 n_end lambda_start lambda_end]); % Set axis limits based on the ranges (0.1 to 1)
% hold off;
% 
% % Write results to an Excel file for n >= 1
% filename_greater_than_1 = 'lambda_n_values_filtered_greater_than_1.xlsx';
% 
% if isfile(filename_greater_than_1)
%     delete(filename_greater_than_1); % Delete the existing file if it exists
% end
% 
% result_table_greater_than_1 = table(lambda_greater_than_1', n_greater_than_1', 'VariableNames', {'Lambda', 'N'});
% writetable(result_table_greater_than_1, filename_greater_than_1);