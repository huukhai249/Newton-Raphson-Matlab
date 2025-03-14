% Constants
k = 0.2; 
beta1_sq = (k^2-1)/(2*log(k));% Given k
beta1 = sqrt(beta1_sq);     % Fixed parameter beta1

% Read data from Excel file (assuming the file is named 'data.xlsx')
data = readmatrix('lambda_n_lt1.xlsx'); % Read the Excel file
beta_values = data(:, 1);       % Column 1: beta values
n_values = data(:, 2);          % Column 2: n values

% Calculate Ks for each pair of beta and n
Ks_values = zeros(size(beta_values)); % Preallocate array for Ks
for i = 1:length(Ks_values)
    n = n_values(i);
    beta = beta_values(i);
    s = 1 / n; % Define s = 1/n
    
    % Compute the first term
    term1 = (1 + k^2 - 2 * beta1^2)/4;   
    part1 = (term1)^(1/(n-1));
    % Compute the second term
    term2_ts = (3+s)*(1-k^2);
    term2_ms_p1 = (1-beta^2)^(1+s);
    term2_ms_p2 = k^(1-s);
    term2_ms_p3 = (beta^2-k^2)^(1+s);

    term2_ms = term2_ms_p1 - term2_ms_p2*term2_ms_p3;
   
    term2 = term2_ts/term2_ms;
    
    part2 = term2^(n/(n-1));

    % Compute Ks
    Ks_values(i) =(1-k)* part1 * part2;
end

% Create the plot
figure;
hold on;

% Use unique beta values for different colors/markers
unique_betas = unique(beta_values);
colors = lines(length(unique_betas)); % Generate distinct colors

plot(n_values, Ks_values);

% Add labels and legend
xlabel('n');
ylabel('K_s');
title('K_s vs n for different β values');
grid on;

% Set axis limits (adjust as needed)
axis([min(n_values)-0.1 max(n_values)+0.1 0 max(Ks_values)*1.1]); % Buffer for visibility
hold off;

ks_mean = mean(Ks_values)

