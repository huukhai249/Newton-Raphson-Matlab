clc; clear;

%% Parameters
k_start = 0.1; % in Xuesi's work, lambda is beta
k_end = 1;
alpha_start = 0.1;
alpha_end = 20.0;
mu = 1; % pa.s
L_in = 0.042; % m - L initial
R_in = 0.006; % m
data = readmatrix('newtonian_pQ_kenic.xlsx');
Q_values = data(:, 2); % Column 2: Q values
DelP_values = data(:, 5); % Column 5: DelP values 
ratioPQ = Q_values ./ -(DelP_values);
mpq = mean(ratioPQ);
fR = (8 * mu * L_in * mpq) / (pi * R_in^4); % sigma_fixed take from find_sigma.m program

%% ----------  START NEWTON-RAPHSON LOOP  --------------%%

k_range = k_start:0.01:k_end; % Range for k
alpha_range = alpha_start:0.01:alpha_end; % Range for alpha (initial guess)

% Tolerances and iteration limit for Newton-Raphson
tolerance = 1e-12; % Adjusted tolerance
max_iterations = 100;

% Initialize arrays to store results
k_values = [];
alpha_values = [];
f_values = [];

fprintf('Iteration\tk\t\tAlpha Init\tAlpha\n');
fprintf('-------------------------------------------------------------\n');

%% Loop over all combinations of alpha and k
for k = k_range
    found = false;
    for alpha_init = alpha_range
        alpha = alpha_init; % Initial guess for alpha
        
        % Newton-Raphson iteration
        for iter = 1:max_iterations
            % Define the function f(alpha) = fx(alpha) - fk(k)
            fk = (1 - k.^2).^2 .* (1 + 1 ./ log(k));
            f_alpha = fR * alpha;
            
            % Function f(k) - f(alpha)
            f = fk - f_alpha;
            df = -fR; % Derivative of f with respect to alpha
            
            % Update alpha
            alpha_new = alpha - f / df;
            
            % Display current guess and function value
            fprintf('%d\t\t%.2f\t%.2f\t%.6f\n', iter, k, alpha_init, alpha_new);
            
            % Check convergence
            if abs(alpha_new - alpha) < tolerance
                % Store the converged values
                k_values = [k_values, k];
                alpha_values = [alpha_values, alpha_new];
                f_values = [f_values, abs(f)]; % Store the absolute value of f(alpha)
                found = true;
                break;
            end
            
            alpha = alpha_new;
            
            % Check if alpha is within bounds (alpha_start to alpha_end)
            if alpha < alpha_start || alpha > alpha_end
                disp('alpha out of bounds, restarting with new initial guess.');
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

%% Create the plot
% Fit a polynomial of degree 1 (linear fit)
% Fit a polynomial of degree 1 (linear fit)
p = polyfit(k_values, alpha_values, 1);

% Generate fitting values using the polynomial
k_fit = linspace(min(k_values), max(k_values), 100);
alpha_fit = polyval(p, k_fit);

% Calculate R-squared
alpha_fit_values = polyval(p, k_values);
SS_res = sum((alpha_values - alpha_fit_values).^2);
SS_tot = sum((alpha_values - mean(alpha_values)).^2);
R_squared = 1 - (SS_res / SS_tot);

% Display R-squared value
fprintf('R-squared: %.4f\n', R_squared);

% Display the polynomial expression
fprintf('Fitting Polynomial: alpha(k) = %.4f*k + %.4f\n', p(1), p(2));

% Plot the original data and the fitting curve
figure;
hold on;
plot(k_values, alpha_values, 'bo', 'DisplayName', 'Newton Raphson Data');
plot(k_fit, alpha_fit, 'r-', 'DisplayName', 'Linear Fit');
xlabel('k');
ylabel('\alpha');
title('The relationship between \alpha and k for Kenic Mixer Geometry');
legend;
grid on;
hold off;