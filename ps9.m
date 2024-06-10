% SPDX-License-Identifier: GPL-3.0-or-later
%
% ps9.m -- Submission for Problem Set 9 in ECE 211-1 session
% Copyright (C) 2024 Jonghyeok Kim <jonghyeok.kim@cooper.edu>

% Clear previous commands
clc;
clear;
close all;

%% Part 1
N = 1e6;
gaussian_data = randn(1, N); % iid N(0, 1)
t_data = trnd(5, 1, N) * sqrt(3 / 5); % As initial variance = 5/3, we need to multiply by sqrt of its reciprocal to have variance of 1
cauchy_data = tan(pi * rand(1, N)) * 0.544;

% Part (a) 
% The fraction of samples where the absolute value is less than 1
disp('1 - (a)');
frac_gaussian = mean(abs(gaussian_data) < 1)
frac_t = mean(abs(t_data) < 1)
frac_cauchy = mean(abs(cauchy_data) < 1)

% Graphing the data values
figure;
sgtitle('1 - (a) Graphs');
subplot(3,1,1);
plot(gaussian_data);
hold on;
xline(-1, 'k--'); % superimpose dashed line at x = -1
xline(1, 'k--'); % superimpose dashed line at x = +1
title('Gaussian Distribution');
xlabel('Sample Index');
ylabel('Amplitude');
hold off;

subplot(3,1,2);
plot(t_data);
hold on;
xline(-1, 'k--');
xline(1, 'k--');
title('Student''s t-Distribution');
xlabel('Sample Index');
ylabel('Amplitude');
hold off;

subplot(3,1,3);
plot(cauchy_data);
hold on;
xline(-1, 'k--');
xline(1, 'k--');
title('Cauchy Distribution');
xlabel('Sample Index');
ylabel('Amplitude');
hold off;

% Part (b)
length_segments = 100000; % length of segments
num_segments = 10; % number of segments
segment_means = zeros(1, num_segments);

for i = 1:num_segments
    segment_start = (i-1) * length_segments + 1; % 1, 100,001, 200,001, ..., 900,001
    segment_end = i * length_segments; % 100,000, 200,000, 300,000, ..., 1,000,000
    segment_means(i) = mean(cauchy_data(segment_start:segment_end));
end

disp('1 - (b)');
disp('Means of each segment for Cauchy distribution:');
disp(segment_means);

%% Part 2
var = 2; % var = (σ_v)^2
disp('2 - (a)');

% 2 - (a) - 1
disp('The model is an ARMA(2,2) model');

% 2 - (a) - 2
disp('It is an innovation filter'); 
% as the filter is transforming white noise v[n] into a signal x[n]

% 2 - (a) - 3
disp('v[n] = x[n]-1.6x[n-1]+0.81x[n-2]-0.4v[n-1]-0.2v[n-2]');

% 2 - (a) - 4
disp('H(z) = (1+0.4z^{-1}+0.2z^{-2})/(1-1.6z^{-1}+0.81z^{-2})');

% 2 - (a) - 5
b = [1, 0.4, 0.2];
a = [1, -1.6, 0.81];
disp('S_x (w) = (σ_v)^2 *|B(w)|^{2}/(|A(w)|^{2})');

% 2 - (a) - 6
[z,p,k] = tf2zp(b,a);

% Checking whether the system is minimum-phase or not
if (b(1) == 0)
    disp('The system is minimum-phase.');
else
    disp('The system is not minimum-phase.');
end

figure;
zplane(z, p);
title('Pole-Zero Plot for 2 - (a) - 6');

% 2 - (b)
disp('');
disp('2 - (b)');
N = 1e4;
v = sqrt(var) * randn(1, N);

% 2 - (b) - (1) Generating x by applying the filter v
x = filter(b, a, v);

% 2 - (b) - (2) Using time-averaging to estimate r_x [m]
max_lag = 6;
r_x = zeros(1, 7);

for m = 0:max_lag
    x1 = x(1:N-m);  % Before being shifted
    x2 = x(m+1:N);  % After being shifted
    
    sum_product = dot(x1, x2);
    r_x(m + 1) = sum_product / (N - m);
end

% 2 - (b) - (3) Displaying the auto-correlation values
% As r_x(m) = r_x(-m), we only need to consider 0=<m=<6 range
full_r_x = [fliplr(r_x(2:end)), r_x];  % A symmetric autocorrelation array
lags = -max_lag:max_lag;
figure;
stem(lags, full_r_x, 'filled');
title('A stem plot of r_x(m) for -6 =< m =< 6');
xlabel('Lags');
ylabel('Autocorrelation value');

% 2 - (b) - (4) Displaying the Toeplitz matrix
R1 = toeplitz(full_r_x) 
disp('Toeplitz matrix R:');
disp(R1);

% 2 - (b) - (5) Compute the eigenvalues
% Using For eigen value λ, q^{H}Aq = λ*(abs(q))^{2}
% Eigen value must be greater than 0 for the matrix to be verified as positive definite.
eigen_values = eig(R1);

isPositiveDefinite = all(eigen_values > 0);
disp(['Matrix R is positive definite: ', num2str(isPositiveDefinite)]);

% 2 - (b) - (6) Arranging the data values x into a Toeplitz matrix
toeplitz_data = toeplitz(x(1:6));

% Estimating the correlation matrix
R2 = toeplitz_data * toeplitz_data'/ N;

% Comparing R1 and R2
norm_dif = norm(R1(1:6,1:6) - R2);
disp('Norm of (R1-R2): ');
disp(norm_dif);

% 2 - (c)
disp('');
disp('2 - (c)');

% 2 - (c) - (1)
[s_est, w] = pwelch(x, hamming(512), 256, 512);
figure;
plot(w, s_est);
title('Power Spectral Density of x');
xlabel('Frequency (rad/sample)');
ylabel('Power/Frequency (units^2/Hz)');

% 2 - (c) - (2) A peak frequency
[~, idx] = max(s_est); % idx = index of max(s_est)
peak_frequency = w(idx);  % Frequency corresponding to the peak
disp(['Peak frequency (ω0) is approximately: ', num2str(peak_frequency), ' rad/sample']);

% 2 - (c) - (3)
pole_angles = abs(angle(p))
disp(['w_0:', num2str(peak_frequency)]);
disp('Thus, the angles are very similar to each other.');

% 2 - (d) AR Modeling
[a, varv] = aryule(x, 4);

% Compare varv to the variance of the innovations signal v (assumed known)
sigma_v2 = 2;  % Example: assuming the variance of v is known and is 2
disp('');
disp(['(σ_v)^2 (=Variance of the innovations signal v): ', num2str(sigma_v2)]);
disp(['varv (=Variance from AR model): ', num2str(varv)]);
disp('Thus, two variances are very similar to each other.');

x0 = filter(1, a, v);
r_x0 = zeros(1, 2*max_lag + 1); % range of m: -6 =< m =< 6
for m = 0:2*max_lag
    x01 = x0(1:N-m);  % Before being shifted
    x02 = x0(m+1:N);  % After being shifted
    
    sum_product0 = dot(x01, x02);
    if m < max_lag
        r_x0(m + 1) = sum_product0 / (N + m);
    else
        r_x0(m + 1) = sum_product0 / (N - m);
    end
end

disp('Autocorrelation of x0:');
disp(r_x0);

% Comparing the first 100 points of x0 and x
figure;
stem(1:100, x(1:100), 'r', 'filled'); 
hold on; 
stem(1:100, x0(1:100), 'b');
title('Comparison of Original and AR Modeled Signals');
xlabel('Sample Index');
ylabel('Amplitude');
legend('Original Signal x', 'AR Modeled Signal x0');
