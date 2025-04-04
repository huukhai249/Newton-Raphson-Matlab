clc;
clear;

%% Process for newtonian fluid
% data = readmatrix('newtonian_pQ_rect.xlsx');
% Q_values = data(:, 3);       % Column 1: beta values
% DelP_values = data(:, 5); 
% % R1 = 0.0025; % m
% % A = pi*R1^2;
% % L = 0.05; 
% % V = A*L;
% %For rectangular cross-section
% b =0.005;
% d =0.005;
% L = 0.05;
% %%%%%%%%%%%
% mu = 1;
% nloop = length(Q_values);
% Kp_values = [];
% for i=1:nloop
% Q = Q_values(i); delP = -DelP_values(i);
% % Kp_values(i) = (delP*Q/V)*((pi^2*R1^6)/(2*Q^2*mu));
% Kp_values(i) = (delP/Q)*(b*d^3)/(36*L*mu);
% end
% Q_values_ml = 10^6*Q_values;
% figure;
% plot(Q_values_ml,Kp_values,'b','LineWidth', 2);
% xlabel('Volume Flowrate,Q [ml/s]');
% ylabel('Kp [-]');
% grid on;
% mean_Kp = mean(Kp_values)
% % Set x and y limits
% ylim([0.4, 1]);
% xlim([-0.3, 2.4]);

%% Process for non-newtonian fluid
data = readmatrix('nonnewtonian_pQ_rect.xlsx');
Q_values = data(:, 4);       % Column 1: beta values
DelP_values = data(:, 5); 
% R1 = 0.0025; % m
% A = pi*R1^2; 
% L = 0.05; 
% For rectangular geometry
b =0.005;
d =0.005;
L = 0.05;
A = b*d;
V = A*L;
rho = 1000;
H = d/6; %H =D/4
kp = 0.7831;
ks = 1.0187;
m =10; n=0.6;

Np_values = zeros(size(Q_values));
Re_values = zeros(size(Q_values));
nloop = length(Q_values);
mu_eff_values = zeros(size(Q_values));
shear_rate_eff_values = zeros(size(Q_values));
% Ks_values = zeros(size(Q_values));

for i = 1:nloop
    Q = Q_values(i);
    Minus_delP = -DelP_values(i);
    u_avg = Q/A;
    shear_rate_app = u_avg/H;

    Np_values(i)=(Minus_delP*Q)/(V*rho*u_avg^2*shear_rate_app);
    Re_values(i) = kp/Np_values(i);

    mu_eff_values(i) = (H^2*Minus_delP*Q)/(kp*V*u_avg^2);
    % shear_rate_eff_values(i) = (mu_eff_values(i)/m)^(1/(n-1)); 
    % Ks_values(i) = shear_rate_eff_values(i)/shear_rate_app;
    shear_rate_eff_values(i) = ks*shear_rate_app;
end


% Plot each value of shear_rate_eff against mu_eff as individual points
figure;

% Q_values_ml = 10^6*Q_values;
% scatter(Q_values_ml,Ks_values,'filled');
% ylim([-10, 10]);
% xlabel('Volume Flowrate,Q [ml/s]');
% ylabel('Ks [-]');
% ks=mean(Ks_values)
% 
% scatter(Re_values,Np_values, 'filled');
% xlabel('Re_e_f_f [-]');
% ylabel('N_p [-]');

scatter(shear_rate_eff_values,mu_eff_values, 'filled');
xlabel('shear-rate_e_f_f [1/s]');
ylabel('\mu_e_f_f [-]');
% 

% 
% Set logarithmic scale
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on;
%%