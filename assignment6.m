% SPDX-License-Identifier: GPL-3.0-or-later
%
% assignment6.m -- Submission for Assignment 6 in ECE 210-B session
% Copyright (C) 2024 Jonghyeok Kim <jonghyeok.kim@cooper.edu>

% Part 2
fs = 192000; % Sampling frequency in Hz
t = (0:fs-1)/(2*fs); % Time vector

frequencies = [-20480, -360, 996, 19840];
amplitudes = [db2mag(14); db2mag(-10); db2mag(0); db2mag(2)];
signals = sum(amplitudes .* sin(2*pi*frequencies.'*t)); % Accumulate signals

noise_level = 10^(-10/20); % White noise
noise = noise_level* randn(size(signals));
signals_with_noise = noise + signals;

% DFT
X = fft(signals_with_noise);

% Plot magnitude on dB scale
figure;
f = (0:length(X)-1)*fs/length(X);
magnitude = 20*log10(abs(X));
plot(f, magnitude)
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on;

% Part 3
k=0.53;
a = [0.76+0.64j; 0.76-0.64j; 0.69+0.71j; 0.69-0.71j; 0.82+0.57j; 0.82-0.57j];
c = [0.57+0.78j; 0.57-0.78j; 0.85+0.48j; 0.85-0.48j; 0.24; 0.64];

figure;
sgtitle('Pole-Zero Plot');

% Zero-pole
subplot(2,2,[1,3]);
zplane(a, c);

% Magnitude Response
subplot(2,2,2);
[zp_b, zp_a] =zp2tf(a,c,k);
w = linspace(0, fs/2, 1e3);
[DE_b, DE_a] = freqz(zp_b, zp_a, w, fs);
mag_DE = 20*log10(abs(DE_b));
phase_DE = unwrap(angle(DE_b));
phaseDeg_DE = phase_DE * 180 / pi;
DE_c = DE_a * fs / (2*pi);
plot(DE_c, mag_DE);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Response');
grid on;

% Phase Response
subplot(2,2,4);
plot(DE_c, phaseDeg_DE);
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase Response');
grid on;

% Part 1
function mag = db2mag(decibel)
    mag = 10^(decibel/20);
end


