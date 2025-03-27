clc; clear;

%% Parameters
beta_start = 0.1; % in Xuesi's work, beta is beta
beta_end = 1;
n_start = 0.1;
n_end = 2.0;
k = 0.25; % Fixed k value (or k value in Xuesi's work)
% k take from find_k.m program
%% ----------  START NEWTON-RAPHSON LOOP  --------------%%

n_range = n_start:0.01:n_end; % Range for n
beta_range = beta_start:0.01:beta_end; % Range for beta (initial guess)

% Tolerances and iteration limit for Newton-Raphson
tolerance = 1e-12; % Adjusted tolerance
max_iterations = 100;
h = 1e-6; % Step size for numerical derivative

% Initialize arrays to store results
n_values = [];
beta_values = [];
f_values = [];

fprintf('Iteration\tk\t\tN\t\tbeta Init\tbeta\n');
fprintf('-------------------------------------------------------------\n');

%% Loop over all combinations of beta and n
for n = n_range
    found = false;
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
            
            % Display current guess and function value
            fprintf('%d\t\t%.6f\t%.2f\t%.2f\t%.6f\n', iter, k, n, beta_init, beta_new);
            
            % Check convergence
            if abs(beta_new - beta) < tolerance
                % Store the converged values
                n_values = [n_values, n];
                beta_values = [beta_values, beta_new];
                f_values = [f_values, abs(f)]; % Store the absolute value of f(beta)
                found = true;
                break;
            end
            
            beta = beta_new;
            
            % Check if beta is within bounds (0.1 to 1)
            if beta < beta_start || beta > beta_end
                disp('beta out of bounds, restarting with new initial guess.');
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

%% Filter the results to get one beta value for each n
filtered_n_values = [];
filtered_beta_values = [];

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
    filtered_beta_values = [filtered_beta_values, beta_values(best_idx)];
end

%% Separate the results into two groups: n < 1 and n >= 1
n_less_than_1 = filtered_n_values(filtered_n_values < 1);
beta_less_than_1 = filtered_beta_values(filtered_n_values < 1);

n_greater_than_1 = filtered_n_values(filtered_n_values >= 1);
beta_greater_than_1 = filtered_beta_values(filtered_n_values >= 1);

%% Create the plot for n < 1
figure; % Create a new figure
hold on; % Allow multiple plots on the same figure

% Plot the results for n < 1
plot(n_less_than_1, beta_less_than_1, 'bo');
xlabel('n');
ylabel('\beta');
title('Relationship between n and \beta for n < 1, fixed k=0.368694, Newtonian fluid');
grid on;

% Ensure the plot is visually appealing
axis([n_start 1 beta_start beta_end]); % Set axis limits based on the ranges (0.1 to 1)
hold off;

% Write results to an Excel file for n < 1
filename_less_than_1 = 'beta_n_lt1.xlsx';

if isfile(filename_less_than_1)
    delete(filename_less_than_1); % Delete the existing file if it exists
end

result_table_less_than_1 = table(beta_less_than_1', n_less_than_1', 'VariableNames', {'beta', 'N'});
writetable(result_table_less_than_1, filename_less_than_1);

beta_mean = mean(beta_less_than_1)


%% Create the plot for n >= 1
% figure; % Create a new figure
% hold on; % Allow multiple plots on the same figure
% 
% % Plot the results for n >= 1
% plot(n_greater_than_1, beta_greater_than_1, 'ro-');
% xlabel('n');
% ylabel('\beta');
% title('Relationship between n and \beta for n >= 1, fixed k=0.368694, Newtonian fluid');
% grid on;
% 
% % Ensure the plot is visually appealing
% axis([1 n_end beta_start beta_end]); % Set axis limits based on the ranges (0.1 to 1)
% hold off;
% 
% % Write results to an Excel file for n >= 1
% filename_greater_than_1 = 'beta_n_values_filtered_greater_than_1.xlsx';
% 
% if isfile(filename_greater_than_1)
%     delete(filename_greater_than_1); % Delete the existing file if it exists
% end
% 
% result_table_greater_than_1 = table(beta_greater_than_1', n_greater_than_1', 'VariableNames', {'beta', 'N'});
% writetable(result_table_greater_than_1, filename_greater_than_1);