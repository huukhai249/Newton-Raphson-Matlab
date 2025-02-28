% Newton-Raphson method for solving f(x) = 0 with additional control
% x is sigma
%% Init Parameters
clc; clear;
data = readmatrix('newtonian_pQ.xlsx');
Q_values = data(:, 3);       % Column 1: beta values
DelP_values = data(:, 6);          % Column 2: n values 
ratio_values = Q_values ./ DelP_values ;
% Constants
L = 200*10e-3;
R = 90/11*10e-3;
ratio = mean(ratio_values);
%% Declare sigma function -- see derivation at ....
A = 8 * L * ratio / (pi * R^4);
f = @(x) ((1 - x.^2).^2) .* (1 + 1 ./ log(x)) - A;

%% Tolerances and iteration limit for Newton-Raphson
tolerance = 1e-16;
max_iterations = 100;
h = 1e-6; % Step size for numerical derivative
control_threshold = 0.0001; % Control threshold for f(x_guess)

% Initial guess
x_guess = 0.3;
guesses_list = zeros(1, max_iterations); % Preallocate for performance
guesses_list(1) = x_guess;

fprintf('Iteration\tGuess\t\tf(Guess)\n');
fprintf('--------------------------------------\n');

%% Newton-Raphson iteration
for i = 1:max_iterations
    x_plus_h = x_guess + h;
    x_minus_h = x_guess - h;

    % Compute f(x + h) and f(x - h)
    f_plus = f(x_plus_h);
    f_minus = f(x_minus_h);

    % Central difference approximation of derivative
    df = (f_plus - f_minus) / (2 * h);

    % Update guess using Newton-Raphson formula
    x_new = x_guess - f(x_guess) / df;

    % Display current guess and function value
    fprintf('%d\t\t%.12f\t%.12f\n', i, x_guess, f(x_guess));

    % Check for convergence
    if abs(x_new - x_guess) < tolerance || abs(f(x_new)) < control_threshold
        guesses_list(i+1) = x_new;
        guesses_list = guesses_list(1:i+1); % Trim unused preallocated space
        break;
    end

    % Update guess
    x_guess = x_new;
    guesses_list(i+1) = x_new;
end

%% Plotting the function and the guesses
t = -1:0.01:3;
y = f(t);

figure;
plot(t, y);
hold on;
plot(guesses_list, f(guesses_list), 'ro');
xlabel('k');
ylabel('F(k)');
title('Newton-Raphson Method to find k');
legend('F(k)', 'Guesses Data');

% Add text box with the equation
dim = [.2 .5 .3 .3]; % Position of the text box
str = 'F(k) = (1-k^2)^2(1+1/ln(k)) - 8\muLQ/\piR^4(-\DeltaP)';
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
grid on;
hold off;