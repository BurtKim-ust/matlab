% SPDX-License-Identifier: GPL-3.0-or-later
%
% assignment1.m -- Submission for Assignment 1 in ECE 210-B session
% Copyright (C) 2024 Jonghyeok Kim <jonghyeok.kim@cooper.edu>

% Part 1
fprintf('Answers for Part 1: ');
% (1)
u = [11, 13, 17];
display(u);

% (2)
v = [-1; -1; -1];
display(v);

% (3)
A = [-u; 2 * u; 7 * u];
display(A);

% (4)
B = [A.', v];
display(B);

% Part 2
fprintf('Answers for Part 2: ');
% (5)
c = exp(1j * pi / 4);
display(c);

% (6)
d = sqrt(1j);
display(d);

% (7)
l = floor(nthroot(8.4108 * 10^6, 2.1));
display(l);

% (8)
k = floor(100 * log(2)) + ceil(exp(7.5858));
display(k);

% Part 3
A = [1, -11, -3
     1, 1, 0
     2, 5, 1];
b = [-37; -1; 10];
x = A \ b;
fprintf('Answer for Part 3: ');
display(x);
