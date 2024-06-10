% SPDX-License-Identifier: GPL-3.0-or-later
%
% ps8.m -- Submission for Problem Set 8 in ECE 211-1 session
% Copyright (C) 2024 Jonghyeok Kim <jonghyeok.kim@cooper.edu>

% Clear previous commands
clc;
clear;
close all;

%% Qustions from Pg 3
disp('ω(θ_1) − ω(θ_2) = 2πk for some integer k');
disp('2π * d / λ_0 * (cos(θ_1) - cos(θ_2)) = 2πk');
disp('This simplifies to d / λ_0 * (cos(θ_1) - cos(θ_2)) = k');
disp('For k == 0, there is no aliasing since θ_1 = θ_2.');
disp('For k =! 0, there is aliasing when value for (d / λ_0) is an integer for θ_1 != θ_2');
disp('The condition for spatial aliasing not to occur is that the value (d / λ_0) must be');
disp('less than 1/2 as this ensures cos(θ_1) - cos(θ_2) does not become an ineger value.');
disp('');
disp('If (d / λ_0) is such that the range of w exceeds [-π, π], there will be values of');
disp('w for which there is no corresponding AOA within the range of [0, π].')
disp('');
disp('(d / λ_0) = 1/2 is the unique value under which no spatial aliasing occurs and there');
disp('is no invisible region.');
disp('The maximum value of w(θ) is π, and the minimum is -π, completely covering the detectable'); 
disp('range without exceeding it. Thus, there exists no invisible region.');
disp('At the same time, this is the unique value where the sensor array can detect signals from');
disp('all possible AOAs without aliasing as it is the unique value that corresponds to the spatial Nyquist rate.');

%% Part I: Signal Generation
M=100;
N=100;
d_lambda = 1/2; % Ratio of d to lambda_0; this is what the signal behavior depends on.
PdB = [0, -2, -4]; % Source powers in dB
Pn_dB = 10; % Noise power in dB
AOAs_1 = [10, 25, 70]; % Scenario 1
AOAs_2 = [10, 12, 70]; % Scenario 2

% Source vector matrix
S_1 = zeros(M, length(AOAs_1)); 
for i = 1:length(AOAs_1)
    theta_1 = AOAs_1(i) * pi / 180; % Convert AOA from degrees to radians
    S_1(:, i) = exp(-1j * 2 * pi * d_lambda * cos(theta_1) * (0:M-1).') / sqrt(M);
end

S_2 = zeros(M, length(AOAs_2)); 
for i = 1:length(AOAs_2)
    theta_2 = AOAs_2(i) * pi / 180;
    S_2(:, i) = exp(-1j * 2 * pi * d_lambda * cos(theta_2) * (0:M-1).') / sqrt(M);
end

% B matrix
variances = 10.^(PdB / 10); % Convert dB to linear scale for variances
B_1 = zeros(length(AOAs_1), N);
for i = 1:length(AOAs_1)
    B_1(i, :) = sqrt(variances(i)) * (randn(1, N) + 1j * randn(1, N)) / sqrt(2); % Matrix p X q = Matrix 1X N in this case
end

B_2 = zeros(length(AOAs_2), N);
for i = 1:length(AOAs_2)
    B_2(i, :) = sqrt(variances(i)) * (randn(1, N) + 1j * randn(1, N)) / sqrt(2);
end

% Noise vector matrix V
noise_var = 10^(Pn_dB / 10);
V = sqrt(noise_var) * (randn(M, N) + 1j * randn(M, N)) / sqrt(2);

% Final composite data matrix A and correlation matrix R
A_1 = S_1 * B_1 + V / sqrt(M);
R_1 = (1/N) * (A_1 * A_1');

A_2 = S_2 * B_2 + V / sqrt(M);
R_2 = (1/N) * (A_2 * A_2');

%% Part II
% Matrices for sval and eigval
[U_1, Sval_1, V_1] = svd(A_1);
[eigvec_1, eigval0_1] = eig(R_1);
[eigval_1, idx_1] = sort(diag(eigval0_1), 'descend');
eigvec_1 = eigvec_1(:, idx_1);

[U_2, Sval_2, V_2] = svd(A_2);
[eigvec_2, eigval0_2] = eig(R_2);
[eigval_2, idx_2] = sort(diag(eigval0_2), 'descend');
eigvec_2 = eigvec_2(:, idx_2);

figure;
sgtitle('Singular values of A and Eigenvalues of R');

subplot(2,2,1);
stem(diag(Sval_1), 'filled');
title('Scenario 1: Singular Values of A');
xlabel('Index');
ylabel('Value');
grid on;

subplot(2,2,2);
stem(eigval_1, 'filled');
title('Scenario 1: Eigenvalues of R');
xlabel('Index');
ylabel('Value');
grid on;

subplot(2,2,3);
stem(diag(Sval_2), 'filled');
title('Scenario 2: Singular Values of A');
xlabel('Index');
ylabel('Value');
grid on;

subplot(2,2,4);
stem(eigval_2, 'filled');
title('Scenario 2: Eigenvalues of R');
xlabel('Index');
ylabel('Value');
grid on;

% The third largest to the fourth largest ratios
svRatio_1 = Sval_1(3,3) / Sval_1(4,4); % Singular Value
evRatio_1 = eigval_1(3) / eigval_1(4); % Eigenvalue
svRatio_2 = Sval_2(3,3) / Sval_2(4,4);
evRatio_2 = eigval_2(3) / eigval_2(4);
disp("Scenario 1: ");
disp(['σ3/σ4 of singular values: ', num2str(svRatio_1)]);
disp(['σ3/σ4 pf eigenvalues: ', num2str(evRatio_1)]);
disp("Scenario 2: ");
disp(['σ3/σ4 of singular values: ', num2str(svRatio_2)]);
disp(['σ3/σ4 pf eigenvalues: ', num2str(evRatio_2)]);

% Projection matrix
L_1 = length(AOAs_1); % Assuming number of sources equals number of AOAs
UL_1 = U_1(:, 1:L_1);
PS_1 = UL_1 * UL_1';
PN_1 = eye(M) - PS_1;

L_2 = length(AOAs_2); % Assuming number of sources equals number of AOAs
UL_2 = U_2(:, 1:L_2);
PS_2 = UL_2 * UL_2';
PN_2 = eye(M) - PS_2;

% MUSIC and MVDR spectrums
theta = 0:0.2:180; % theta range: 0 degree - 180 degree
SMUSIC_1 = zeros(size(theta));
SMVDR_1 = zeros(size(theta));
SMUSIC_2 = zeros(size(theta));
SMVDR_2 = zeros(size(theta));
Rinv_1 = inv(R_1);
Rinv_2 = inv(R_2);

for i = 1:length(theta)
    AOA_1 = theta(i) * pi / 180;
    s_1 = exp(-1j * 2 * pi * d_lambda * cos(AOA_1) * (0:M-1).') / sqrt(M);
    SMUSIC_1(i) = 1 / (s_1' * PN_1 * s_1);
    SMVDR_1(i) = 1 / (s_1' * Rinv_1 * s_1);
end

for i = 1:length(theta)
    AOA_2 = theta(i) * pi / 180;
    s_2 = exp(-1j * 2 * pi * d_lambda * cos(AOA_2) * (0:M-1).') / sqrt(M);
    SMUSIC_2(i) = 1 / (s_2' * PN_2 * s_2);
    SMVDR_2(i) = 1 / (s_2' * Rinv_2 * s_2);
end

figure;
sgtitle('S\_MUSIC and S\_MVDR');

% Scenario 1: S_MUSIC
subplot(2,2,1);
title('Scenario 1: S\_MUSIC(θ)');
plot(theta, real(SMUSIC_1));
xlabel('\theta (degrees)');
ylabel('S\_MUSIC');
grid on;

% Scenario 1: S_MVDR
subplot(2,2,2);
title('Scenario 1: S\_MVDR(θ)');
plot(theta, real(SMVDR_1));
xlabel('\theta (degrees)');
ylabel('S\_MVDR');
grid on;

% Scenario 2: S_MUSIC
subplot(2,2,3);
title('Scenario 2: S\_MUSIC(θ)');
plot(theta, real(SMUSIC_2));
xlabel('\theta (degrees)');
ylabel('S\_MUSIC');
grid on;

% Scenario 2: S_MVDR
subplot(2,2,4);
title('Scenario 2: S\_MVDR(θ)');
plot(theta, real(SMVDR_2));
xlabel('\theta (degrees)');
ylabel('S\_MVDR');
grid on;

% For the comments
disp(' '); 
disp('Comments for Part II: ');
disp('Yes, I can see clear peaks');
disp('Yes, I agree that MUSIC is better than MVDR because the graphs generated by');
disp('MUSIC have more clear peaks than the graphs generated by MVDR.');
disp('It seems like MUSIC has a greater ability in sensing the noise and ');
disp('providing useful data(peaks) in more visible way.');
%% Part III
SHS_1 = S_1' * S_1;
abs_SHS_1 = abs(SHS_1);
disp(' '); 
disp('Part III: ');
disp('Scenario 1: Matrix of absolute values of S^{H}S: ');
disp(abs_SHS_1);

SHS_2 = S_2' * S_2;
abs_SHS_2 = abs(SHS_2);
disp('Scenario 2: Matrix of absolute values of S^{H}S: ');
disp(abs_SHS_2);

% Explanation to the questions
disp(' '); 
disp('Comments for Part III: ');
% The SHS matrix reveals the correlation between steering vectors.
disp('Diagonal elements are 1s because of the self-correlation.');
disp('Off-diagonal elements show the correlation between different source vectors.');
disp('Based on the results from scenarioes 1 and 2, the values closer to 1 indicate');
disp('higher similarity with other values.');
disp('For example, 70 degree is more differentiable from the second term in ');
disp('scenario 2(second term: 12) than scenario 1(second term: 25).');
disp('Thus, it led to the lower correlation value in scenario 2(=0.0068)');
disp('than the value in scenario 1(=0.0081).');
disp('Therefore, these matrices provide how correlate their values are to other values.');