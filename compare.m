clc;
clear;

% Read data from the Excel file
data = readtable('compare.xlsx');

% Extract the relevant columns from the table
Np_rec = data{:,9};
Re_rec = data{:,10};
mu_rec = data{:,11};
sr_rec = data{:,12};

Np_kenic = data{:,5};
Re_kenic = data{:,6};
mu_kenic = data{:,7};
sr_kenic = data{:,8};

Np_A = data{:,13};
Re_A = data{:,14};
mu_A = data{:,15};
sr_A = data{:,16};


f = @(y)(10*y.^(-0.4));
y = linspace(0.001,100, 50);
% Create a scatter plot for the data
% figure;
% scatter(Re_A, Np_A, 'filled');
% hold on;
% scatter(Re_rec, Np_rec, 'filled');
% % Add labels to the axes
% xlabel('Re_{eff} [-]');
% ylabel('N_p [-]');
% %Add a legend to distinguish between the two data sets
% legend('Annulus Geometry', 'Original Geometry');



% scatter(sr_A, mu_A, 'filled');
% hold on;
% scatter(sr_rec, mu_rec, 'filled');
% scatter(sr_kenic, mu_kenic, 'filled');
% 
% plot(y,f(y));
% % Add labels to the axes
% xlabel('Shear-rate_{eff} [1/s]');
% ylabel('\mu_{eff} [-]');
% 
% % Add a legend to distinguish between the two data sets
% legend('Annulus Geometry', 'Rectangular Geometry','Kenic Geometry','Theorical curve');
% 
% set(gca, 'XScale', 'log');
% % Set y-axis to logarithmic scale
% set(gca, 'YScale', 'log');
% grid on;
%% Compute error 
% From Rectangular Jang frame work 
np = length(mu_rec);
mu_error = [];
mu_theo = [];

for i =1:np
    mu_rect = mu_rec(i);
    sr_rect = sr_rec(i);
    mu_theo_val = 10*sr_rect^(-0.4);
    mu_error(i) = 100*(mu_rect - mu_theo_val)/mu_theo_val;
    mu_theo(i)= mu_theo_val;
end

scatter(mu_theo,mu_error);
% set(gca, 'XScale', 'log');
% Set y-axis to logarithmic scale
% set(gca, 'YScale', 'log');
grid on;
axis([min(mu_theo)-1 max(mu_theo)+1 min(mu_error)-0.5 max(mu_error)+0.5]); % Buffer for visibility
xlabel('Shear-rate theory rectangular [1/s]');
ylabel('Error [%]');
% title('Error between the theoretical and simulation solutions for Annulus geometry');

hold off;
