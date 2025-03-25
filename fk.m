clc;clear;

% Định nghĩa khoảng giá trị của k
k = linspace(0.01, 0.999, 1000); % Tránh k = 0 vì ln(0) không xác định

% Tính toán giá trị của hàm f(k)
f_k = (1 - k.^2).^2 .* (1 + 1 ./ log(k));

% Vẽ đồ thị
figure;
plot(k, f_k, 'b-', 'LineWidth', 2);
hold on;

% Vẽ đường y = 0
yline(0.8, 'red--', 'LineWidth', 1.5);
yline(0.5, 'black--', 'LineWidth', 1.5);
yline(0.2, 'green--', 'LineWidth', 1.5);

% Thêm nhãn và tiêu đề
xlabel('k');
ylabel('f(k)');
grid on;
legend('f(k)');

% Thiết lập trục
xlim([0 1.2]); % Giới hạn trục x từ 0 đến 2
ylim([-0.6 1.1]); % Giới hạn trục y
hold off;