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
figure;
scatter(Re_A, Np_A, 'filled');
hold on;
scatter(Re_rec, Np_rec, 'filled');
% Add labels to the axes
xlabel('Re_{eff} [-]');
ylabel('N_p [-]');
%Add a legend to distinguish between the two data sets
legend('Annulus Geometry', 'Original Geometry');



% scatter(sr_A, mu_A, 'filled');
% hold on;
% scatter(sr_rec, mu_rec, 'filled');
% scatter(sr_kenic, mu_kenic, 'filled');
% 
% plot(y,f(y));
% % Add labels to the axes
% xlabel('Shear-rate_{eff} [1/s]');
% ylabel('\mu_{eff} [-]');

% Add a legend to distinguish between the two data sets
% legend('Annulus Geometry', 'Rectangular Geometry','Kenic Geometry','Theorical curve: \mu=10\dot{\gamma}');
% 
% hold off;
set(gca, 'XScale', 'log');
% Set y-axis to logarithmic scale
set(gca, 'YScale', 'log');
grid on;