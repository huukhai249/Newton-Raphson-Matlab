clc; clear;

%% Read simulation data

data = readmatrix('newtonian_pQ_kenic.xlsx');
Q_values = data(:, 2);       % Column 1: beta values
DelP_values = data(:, 5); 

mu = 1; %user defind in newtonian simulation 
A_avg = 450e-6; % m^2
V = 9e-5; % m^3
H = 180/11*10e-4; % m
rho = 1000;
kp = 38.3338;
n = 0.7; % power-law index
m = 10; % powerlaw factor

% Np_values = zeros(size(Q_values));
% kp_values = zeros(size(Q_values));
% Re_values = zeros(size(Q_values));
Np_values = [];
kp_values = [];
Re_values = [];
shear_rate_eff_values = []; % Preallocate array for shear_rate_eff
mu_eff_values = [];
nstep = size(Q_values);

%% LOOP FOR NEWTONIAN FLUIDS
% for i = 1:nstep
% Q = Q_values(i);
% delP = -DelP_values(i); 
% u_avg = Q/A_avg;
% y_app = u_avg/H;
% Re_values(i) = rho*u_avg*H/mu;
% kp_values(i) = (delP*Q)/(V*mu*y_app^2);
% Np_values(i) = (delP*Q)/(V*rho*u_avg^2*y_app);
% end
% 
% kp = mean(kp_values)
%% LOOP FOR NON-NEWTONIAN FLUIDS
for i = 1:nstep
Q = Q_values(i);
delP = -DelP_values(i); 
u_avg = Q/A_avg;
y_app = u_avg/H;
Np_values(i) = (delP*Q)/(V*rho*u_avg^2*y_app);
mu_eff_values(i) = (H^2*delP*Q)/(kp*V*u_avg^2);
Re_values(i) = (rho*u_avg*H)/mu_eff_values(i);
shear_rate_eff_values(i) = (mu_eff_values(i)/m)^(1/(n-1));
end
%% plot data
scatter(Re_values, Np_values, 'filled');
% axis([min(shear_rate_eff_values)-0.1 max(shear_rate_eff_values)+0.1 0 max(mu_eff_values)*1.1]); % Buffer for visibility

xlabel('Re_e_f_f [-]');
ylabel('N_p [-]');
set(gca, 'XScale', 'log');
% Set y-axis to logarithmic scale
set(gca, 'YScale', 'log');
grid on;