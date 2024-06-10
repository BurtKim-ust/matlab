% SPDX-License-Identifier: GPL-3.0-or-later
%
% ps5.m -- Submission for Problem Set 5 Q1 in ECE 211-1 session
% Copyright (C) 2024 Jonghyeok Kim <jonghyeok.kim@cooper.edu>

% Clear previous commands
clc;
clear;
close all;

% Given Values
fs = 10e6;
fn = fs /2;
Rp = 2; % [dB] value
Wp = [1.5e6, 2e6];
Ws = [1.4e6, 2.2e6];
Wpd = [1.5e6, 2e6] / fn; % Wpd = fdpass = the pass frequency specifications normalized with respect to the Nyquist bandwidth 
Wsd = [1.4e6, 2.2e6] / fn; % Wsd = fdstop = the stop frequency specifications normalized with respect to the Nyquist bandwidth
Rs = 40; % [dB] value

% Parts (a) and (b)
% Analog IIR form of elliptic
[analog_elliptic_n, analog_elliptic_wn] = ellipord(Wp, Ws, Rp, Rs, 's');
[analog_elliptic_z, analog_elliptic_p, analog_elliptic_k] = ellip(analog_elliptic_n, Rp, Rs, analog_elliptic_wn, 's');
fprintf('Analog Elliptic Order: %d\n', 2 * analog_elliptic_n);

% Digital IIR form of elliptic
[digital_elliptic_n, digital_elliptic_wn] = ellipord(Wpd, Wsd, Rp, Rs);
[digital_elliptic_z, digital_elliptic_p, digital_elliptic_k] = ellip(digital_elliptic_n, Rp, Rs, digital_elliptic_wn);
fprintf('Digital Elliptic Order: %d\n', 2 * digital_elliptic_n);

% Analog IIR form of Chebyshev Type I
[analog_chebyshev_n, analog_chebyshev_wn] = cheb1ord(Wp, Ws, Rp, Rs, 's');
[analog_chebyshev_z, analog_chebyshev_p, analog_chebyshev_k] = cheby1(analog_chebyshev_n, Rp, analog_chebyshev_wn, 's');
fprintf('Analog Chebyshev Type I Order: %d\n', 2 * analog_chebyshev_n);

% Digital IIR form of Chebyshev Type I
[digital_chebyshev_n, digital_chebyshev_wn] = cheb1ord(Wpd, Wsd, Rp, Rs);
[digital_chebyshev_z, digital_chebyshev_p, digital_chebyshev_k] = cheby1(digital_chebyshev_n, Rp, digital_chebyshev_wn);
fprintf('Digital Chebyshev Type I Order: %d\n', 2 * digital_chebyshev_n);

% Parts (c) and (d)
% Analog IIR form of elliptic
figure;
sgtitle('Analog Elliptic Filter Poles and Zeros');

subplot(10, 1, 1); % To avoid overlapping of sgtitle and title
axis off;
subplot(10, 1, [2,10])

zplane(analog_elliptic_z, analog_elliptic_p);

figure;
sgtitle('Analog Elliptic Filter Frequency Response');
subplot(2,1,1); % Magnitude graph
[b, a] = zp2tf(analog_elliptic_z, analog_elliptic_p, analog_elliptic_k);
[AE_b, AE_a] = freqs(b, a, 1e3);
m_AE = 20 * log10(abs(AE_b));
phase_AE = unwrap(angle(AE_b));
phaseDeg_AE = phase_AE * 180 / pi;
plot(AE_a, m_AE);
grid on;
xlim([0, fn]); % limits from DC to the Nyquist bandwidth
ylim([-100, 10]);
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');

subplot(2,1,2); % Phase graph
plot(AE_a, phaseDeg_AE);
grid on;
xlim([0, fn]); % limits from DC to the Nyquist bandwidth
xlabel('Frequency (MHz)');
ylabel('Phase [Degrees]')

% Digital IIR form of elliptic
figure;
sgtitle('Digital Elliptic Filter Poles and Zeros');

subplot(10, 1, 1); % To avoid overlapping of sgtitle and title
axis off;
subplot(10, 1, [2,10])

zplane(digital_elliptic_z, digital_elliptic_p);

figure;
sgtitle('Digital Elliptic Filter Frequency Response');
subplot(2,1,1); % Magnitude graph
[b, a] = zp2tf(digital_elliptic_z, digital_elliptic_p, digital_elliptic_k);
[DE_b, DE_a] = freqz(b, a, linspace(0, pi, 1e3));
m_DE = 20 * log10(abs(DE_b));
phase_DE = unwrap(angle(DE_b));
phaseDeg_DE = phase_DE * 180 / pi;
DE_c = DE_a * fn / pi;
plot(DE_c, m_DE);
grid on;
xlim([0, fn]); % limits from DC to the Nyquist bandwidth
ylim([-100, 10]);
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');

subplot(2,1,2); % Phase graph
plot(DE_c, phaseDeg_DE);
grid on;
xlim([0, fn]); % limits from DC to the Nyquist bandwidth
xlabel('Frequency (MHz)');
ylabel('Phase [Degrees]')

% Analog IIR form of Chebyshev Type I
figure;
sgtitle('Analog Chebyshev Type I Filter Poles and Zeros');

subplot(10, 1, 1); % To avoid overlapping of sgtitle and title
axis off;
subplot(10, 1, [2,10])

zplane(analog_chebyshev_z, analog_chebyshev_p);

figure;
sgtitle('Analog Chebyshev Type I Filter Frequency Response');
subplot(2,1,1); % Magnitude graph
[b, a] = zp2tf(analog_chebyshev_z, analog_chebyshev_p, analog_chebyshev_k);
[AC_b, AC_a] = freqs(b, a, 1e3);
m_AC = 20 * log10(abs(AC_b));
phase_AC = unwrap(angle(AC_b));
phaseDeg_AC = phase_AC * 180 / pi;
plot(AC_a, m_AC);
grid on;
xlim([0, fn]); % limits from DC to the Nyquist bandwidth
ylim([-400, 10])
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');

subplot(2,1,2); % Phase graph
plot(AC_a, phaseDeg_AC);
grid on;
xlim([0, fn]); % limits from DC to the Nyquist bandwidth
xlabel('Frequency (MHz)');
ylabel('Phase [Degrees]')

% Digital IIR form of Chebyshev Type I
figure;
sgtitle('Digital Chebyshev Type I Filter Poles and Zeros');

subplot(10, 1, 1); % To avoid overlapping of sgtitle and title
axis off;
subplot(10, 1, [2,10])

zplane(digital_chebyshev_z, digital_chebyshev_p);

figure;
sgtitle('Digital Chebyshev Type I Filter Frequency Response');
subplot(2,1,1); % Magnitude graph
[b, a] = e(digital_chebyshev_z, digital_chebyshev_p, digital_chebyshev_k);
[DC_b, DC_a] = freqz(b, a, linspace(0, pi, 1e3));
m_DC = 20 * log10(abs(DC_b));
phase_DC = unwrap(angle(DC_b));
phaseDeg_DC = phase_DC * 180 / pi;
DC_c = DE_a * fn / pi;
plot(DC_c, m_DC);
grid on;
xlim([0, fn]); % limits from DC to the Nyquist bandwidth
ylim([-400, 10]);
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');

subplot(2,1,2); % Phase graph
plot(DC_c, phaseDeg_DC);
grid on;
xlim([0, fn]); % limits from DC to the Nyquist bandwidth
xlabel('Frequency (MHz)');
ylabel('Phase [Degrees]')
