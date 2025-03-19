clc;clear;
k=0.3;
data = readmatrix('newtonian_pQ_kenic.xlsx');
Q_values = data(:, 2);       % Column 1: beta values
DelP_values = data(:, 5); 
Kp_values = [];

R = 0.006;
L = 0.042;
V = pi*R^2*L;
mu = 1;


for i = 1:length(Q_values)
minus_DelP = -DelP_values(i);
Q = Q_values(i);

p1 = minus_DelP*Q/V;
p2 = (pi^2*R^6)/(2*Q^2*mu);

Kp_values(i) =  p1*p2;

end

plot(Q_values, Kp_values,'blue','LineWidth',2);
% Set y-axis to logarithmic scale
xlabel('k');
ylabel('K_p');
% title('K_p vary by k');
grid on;

% Set axis limits (adjust as needed)
% axis([min(k_values)-0.1 max(k_values)+0.1 1 25]); % Buffer for visibility
hold off;


mean(Kp_values)
grid on;