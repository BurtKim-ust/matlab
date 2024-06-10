% SPDX-License-Identifier: GPL-3.0-or-later
%
% assignment5.m -- Submission for Assignment 5 in ECE 210-B session
% Copyright (C) 2024 Jonghyeok Kim <jonghyeok.kim@cooper.edu>

clc;
clear all;

% Part 1
syms n t;
n = 0:50;
a_n = 2 .* n + 1;
t = linspace(-pi, pi, 1000);
s = sin(a_n' .* t) ./ a_n';
figure;
hold on;
plot(t, s);
plot(t, sum(s));
title('Part 1')
xlabel('t');
ylabel('Value of s');
xlim([-pi, pi]);
ylim([-1.5, 1.5]);
xticks([-pi, -pi/2, 0, pi/2, pi]);
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});

% Part 2
figure;
hold on;
subplot(1,2,1);
plot(t, s);
title('Components')

subplot(1,2,2);
plot(t, sum(s));
title('Approximation')
sgtitle('Part 2')

% Part 3
syms x y;
x = linspace(-2*pi, 2*pi, 100);
y = linspace(-2*pi, 2*pi, 100);
[X,Y] = meshgrid(x,y);
Z = X .*sin(X)- Y .*cos(Y);
figure;
title('Part 3')
surf(X,Y,Z);
c = parula;
colormap(c);
xticks([-2*pi, -pi, 0, pi, 2*pi]);
xticklabels({'-\2pi', '-\pi', '0', '\pi', '\2pi'});