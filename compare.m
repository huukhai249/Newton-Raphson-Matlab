clc;
clear;

% Read data from the Excel file
data = readtable('compare.xlsx');

% Extract the relevant columns from the table
Np_A = data{:,1};
Re_A = data{:,2};
mu_A = data{:,3};
sr_A = data{:,4};
Np_J = data{:,5};
Re_J = data{:,6};
mu_J = data{:,7};
sr_J = data{:,8};
mu_lw = data{:,9};
sr_lw = data{:,10};


f = @(y)(10*y.^(-0.4));
y = linspace(0.01,100, 50);
% Create a scatter plot for the data
figure;
% scatter(Re_A, Np_A, 'filled');
% hold on;
% scatter(Re_J, Np_J, 'filled');
% % Add labels to the axes
% xlabel('Re_{eff} [-]');
% ylabel('N_p [-]');
% %Add a legend to distinguish between the two data sets
% legend('Annulus Geometry', 'Original Geometry');



scatter(sr_A, mu_A, 'filled');
hold on;
scatter(sr_J, mu_J, 'filled');
scatter(sr_lw, mu_lw, 'filled');

plot(y,f(y));
% Add labels to the axes
xlabel('Shear-rate_{eff} [1/s]');
ylabel('\mu_{eff} [-]');

% Add a legend to distinguish between the two data sets
legend('Annulus Geometry', 'Kenics Geometry','Ansymetric Geometry','Theorical curve');

hold off;
set(gca, 'XScale', 'log');
% Set y-axis to logarithmic scale
set(gca, 'YScale', 'log');
grid on;