clc; clear;

% Read data from Excel file (assuming the file is named 'non_newtonian_pQ.xlsx')
data = readmatrix('non_newtonian_pQ.xlsx'); % Read the Excel file
Q_values = data(:, 2);       % Column 2: Q values
DelP_values = data(:, 5);    % Column 5: DelP values

A_avg = 0.5 * (10 * 30 + 10 * 60) * 10e-6;
P_wet = 2 * (10 + 45) * 10e-3;
L = 200 * 10e-3;
R = 2 * A_avg / P_wet;
k = 0.368694; % Taken from find_sigma.m file
beta1_sq = (k^2 - 1) / (2 * log(k)); 
kp = (8 * (1 - k)^2) / (1 + k^2 - 2 * beta1_sq);
ks = 4.5284;
V = pi * R^2 * L * (1 - k^2);
rho = 1000;
H = R * (1 - k);

shear_rate_eff_values = zeros(size(Q_values)); % Preallocate array for shear_rate_eff
mu_eff_values = zeros(size(Q_values)); % Preallocate array for mu_eff
U_avg_values = zeros(size(Q_values));
Np_values = zeros(size(Q_values));
Re_values = zeros(size(Q_values));

for i = 1:length(Q_values)
    Q = Q_values(i);
    delP = -DelP_values(i);
    shear_rate_app = Q / (pi * R^3 * (1 - k)^2 * (1 + k));
    shear_rate_eff_values(i) = ks * shear_rate_app;
    mu_eff_values(i) = (delP * Q) / (V * kp * shear_rate_app^2);

    u_avg = Q / A_avg;
    Np_values(i) = delP * Q / (rho * u_avg^2 * shear_rate_app);
    Re_values(i) = rho * u_avg * H / (mu_eff_values(i));
end

% Plot each value of shear_rate_eff against mu_eff as individual points
figure;
scatter(shear_rate_eff_values, mu_eff_values, 'filled');
xlabel('Effective Shear rate [1/s]');
ylabel('Effective Viscosity [Pa.s]');
title('Relationship between \mu effective and Shear rate');
% Set x-axis to logarithmic scale
set(gca, 'XScale', 'log');
% Set y-axis to logarithmic scale
set(gca, 'YScale', 'log');
grid on;