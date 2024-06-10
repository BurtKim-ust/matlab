% Define the target correction rate
target_correction_rate = 0.99;

% Define a range of sigma values to evaluate
sigma_values = linspace(0.1, 1.5, 100);

% Preallocate the probabilities array
probabilities = zeros(size(sigma_values));

% Loop over each sigma value to calculate the cumulative error probability
for i = 1:length(sigma_values)
    sigma = sigma_values(i);
    P_error = 0.5 * (erfc(-0.5 / (sigma * sqrt(2))) - erfc((0.5 - 1) / (sigma * sqrt(2))));
    probabilities(i) = sum(binopdf(0:2, 10, P_error));
end

% Plot the probabilities against sigma values
figure;
plot(sigma_values, probabilities, 'b-', 'LineWidth', 2);
hold on;
line([min(sigma_values), max(sigma_values)], [target_correction_rate, target_correction_rate], ...
     'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
title('Cumulative Error Probability vs. Noise Standard Deviation');
xlabel('Noise Standard Deviation (\sigma)');
ylabel('Cumulative Error Probability');
legend('Cumulative Error Probability', 'Target 99% Correction Rate', 'Location', 'SouthEast');
grid on;
