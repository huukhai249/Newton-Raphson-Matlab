clc; clear;

% Read data from Excel file (assuming the file is named 'non_newtonian_pQ.xlsx')
data = readmatrix('nonnewtonian_pQ_annulus_rect.xlsx'); % Read the Excel file
Q_values = data(:, 3);       % Column 2: Q values
DelP_values = data(:, 5);    % Column 5: DelP values
%% for Annulus geometry
k = 0.25;
beta1_sq = (k^2 - 1) / (2 * log(k)); 
kp = (8 * (1 - k)^2) / (1 + k^2 - 2 * beta1_sq);
ks = 4.5162;
L_in = 0.05; % m - L initial
alpha = 0.342625;
L = alpha*L_in;
R= 0.003; % m
A_avg = pi*R^2*(1-k^2);
V = pi * R^2 * L * (1 - k^2);
rho = 1000;
H = R * (1 - k);
%%
shear_rate_eff_values = zeros(size(Q_values)); % Preallocate array for shear_rate_eff
mu_eff_values = zeros(size(Q_values)); % Preallocate array for mu_eff
U_avg_values = zeros(size(Q_values));
Np_values = zeros(size(Q_values));
Re_values = zeros(size(Q_values));

for i = 1:length(Q_values)
    Q = Q_values(i);
    minus_delP = -DelP_values(i);
    shear_rate_app = Q / (pi * R^3 * (1 - k)^2 * (1 + k));
    shear_rate_eff_values(i) = ks * shear_rate_app;
    mu_eff_values(i) = (minus_delP * Q) / (V * kp * shear_rate_app^2);

    u_avg = Q / A_avg;
    Np_values(i) = (minus_delP * Q) / (rho * u_avg^2 * shear_rate_app);
    Re_values(i) = (rho * u_avg * H) / (mu_eff_values(i));
end

% Plot each value of shear_rate_eff against mu_eff as individual points
figure;
% scatter(Re_values,Np_values, 'filled');
% xlabel('Re_e_f_f [-]');
% ylabel('N_p [-]');

scatter(shear_rate_eff_values,mu_eff_values,'filled');
xlabel('shear-rate_e_f_f [1/s]');
ylabel('\mu_e_f_f [-]');
% ratio = mu_eff_values ./ shear_rate_eff_values;
% scatter(shear_rate_eff_values,ratio, 'filled');
% xlabel('shear-rate_e_f_f [1/s]');
% ylabel('ratio [-]');

% scatter(shear_rate_eff_values,mu_eff_values,'filled');
% xlabel('shear-rate_e_f_f [1/s]');
% ylabel('\mu_e_f_f [-]');

% title('Relationship between Re_e_f_f and N_p');
% Set x-axis to logarithmic scale
set(gca, 'XScale', 'log');
% Set y-axis to logarithmic scale
set(gca, 'YScale', 'log');
grid on;