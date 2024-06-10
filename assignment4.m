% SPDX-License-Identifier: GPL-3.0-or-later
%
% assignment4.m -- Submission for Assignment 4 in ECE 210-B session
% Copyright (C) 2024 Jonghyeok Kim <jonghyeok.kim@cooper.edu>

% Part 1
% An anonymous function that performs the standard inner product over ℂ 
ip = @(x, y) sum(x .* conj(y));
% An anonymous function for its associated L^2 norm.
ip_norm = @(x) sqrt(ip(x, x));

% Part 4
% Calculate orthonormal vectors from the set S
S = {complex([1; 1j; 2-1j; -1]), complex([2+3j;3j;1-1j;2j]), complex([-1+7j; 6+10j;11-4j;3+4j])};
U = cell(1, length(S));
% Run gram_schmidt for each of the elements in the set S.
for i = 1:length(S)
    V = gram_schmidt(S{i}, ip, ip_norm);
    U{i} = V;
end

% Part 5
% Initialize scalar orthogonal
orthogonal = true;
% Check orthogonality
for x = 1:length(U)
    [m, n] = size(U{x});
    for i = 1:n
        for j = i+1:n
            if ~isorthogonal(U{x}(:, i), U{x}(:, j), ip)
                orthogonal = false;
                break;
            end
        end
        if ~orthogonal
            break;
        end
    end
end

% Part 2
% A gram_schmidt function
function B = gram_schmidt(A, ip, ip_norm)
    [m, n] = size(A);
    % B = matrix for phi_n (orthonormal set)
    B = zeros(m, n);  
    
    % Gram-Schmidt process
    for j = 1:n
        % v = copy of the jth column vector of A
        v = A(:, j);  
        for i = 1:j-1
            v = v - ip(B(:, i), A(:, j)) / ip_norm(B(:, i))^2 * B(:, i);
        end
        % Adding new component into matrix B 
        B(:, j) = v / ip_norm(v);
    end
end

% Part 3
% An isorthogonal function 
function result = isorthogonal(vector_A, vector_B, ip)
    ip_value = ip(vector_A, vector_B);
    if abs(ip_value) < eps
        result = true;
    else
        result = false;
    end
end
