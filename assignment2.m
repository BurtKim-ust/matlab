% SPDX-License-Identifier: GPL-3.0-or-later
%
% assignment2.m -- Submission for Assignment 2 in ECE 210-B session
% Copyright (C) 2024 Jonghyeok Kim <jonghyeok.kim@cooper.edu>

% Part 1
fprintf("Part 1");
u = -4:2:4
v = 0:pi/4:pi

% Part 2
fprintf("Part 2");
f = prod([1:1:10])

% Part 3-(a)
fprintf("Part 3-(a)");
A = zeros([2 4]);
A(1:1) = 1;
A(2, 3) = 1

% Part 3-(b)
fprintf("Part 3-(b)");
B = reshape([1:2:15,2:2:16], [4,4])

% Part 4 Question
fprintf("Part 4");
a_n = 2 * (0:1:50) + 1;
t = linspace(-pi, pi, 1000);
s = sum(sin(a_n' .* t) ./ a_n');
plot(t, s);


