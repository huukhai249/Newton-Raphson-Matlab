clc; clear;

% Read data from Excel file (assuming the file is named 'data.xlsx')
data = readmatrix('newtonian_pQ.xlsx'); % Read the Excel file
Q_values = data(:, 3);       % Column 1: beta values
DelP_values = data(:, 6);          % Column 2: n values

A_avg = 0.5*(10*30+10*60);
P_wet = 2*(10+45);
L = 200;
R = 2*A_avg/P_wet;
k = 0.367879;
beta1 = 0.657519; 
kp = (8*(1-k)^2)/(1+k^2-2*beta1^2);
ks = 4.5058;
V = pi*R^2*L*(1-k^2);

shear_rate_eff_values = zeros(size(Q_values)); % Preallocate array for Ks
mu_eff_values = zeros(size(Q_values)); % Preallocate array for Ks

for i = 1:length(Q_values)
    Q = Q_values(i);
    delP = -DelP_values(i);
    shear_rate_app = Q/(pi*R^3*(1-k)^2*(1+k));
    shear_rate_eff_values(i) = ks*shear_rate_app;
    mu_eff_values(i) = (delP*Q)/(V*kp*shear_rate_app^2);
end

plot(shear_rate_eff_values,mu_eff_values)