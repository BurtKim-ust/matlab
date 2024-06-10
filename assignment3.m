% SPDX-License-Identifier: GPL-3.0-or-later
%
% assignment3.m -- Submission for Assignment 3 in ECE 210-B session
% Copyright (C) 2024 Jonghyeok Kim <jonghyeok.kim@cooper.edu>

ITERATIONS = 1e6; 
CREWMATES = 6;
ROUNDS = 12;

CREWMATE_SIDES = 4;
IMPOSTER_ROLLS = 2;
IMPOSTER_SIDES = 2;
rng(0x73757300);

% Part 1
crewmates = randi(CREWMATE_SIDES, CREWMATES, ITERATIONS);
sus = sum(randi(IMPOSTER_SIDES, IMPOSTER_ROLLS, ITERATIONS));
targets = randi(CREWMATES, ROUNDS, ITERATIONS);

% Part 2
% Initially setting values in the kills matrix as 0
kills = false(CREWMATES, ITERATIONS);
% Making column for linear Indices
column = repmat(1:ITERATIONS, ROUNDS, 1);
% Use sub2ind to find linear indices representing the positions where kills occur
% Example of using sub2ind: "ind = sub2ind(sz,row,col)"
linearIndices= sub2ind(size(kills), targets, column);
kills(linearIndices) = true;

% Part 3
% Find those who are dead and use that information to find who survivors are
dead = crewmates < sus & kills;
survivors = 1 - dead;

% Part 4
loss_rate = mean(sum(survivors)<2); % Mean of all crewmates surviving across iterations
fprintf('Loss Rate: %.4f\n', loss_rate);