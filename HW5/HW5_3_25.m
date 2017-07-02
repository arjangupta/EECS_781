% Author: Arjan Gupta
% Description: HW#5, 3.25

clear; close all; clc;

% Initialize data-set size
n = 9;
% Assign interpolation points and values at those points
f = [20, 40, 100, 400, 1250, 4000, 16000, 40000, 80000];
A = [.008, .03, .151, .592, 1.477, 9.618, 122.278, 429.310, 850.536];
f_log = log10(f);
A_log = log10(A);

% Initialize extra data points
f_extra = [63, 200, 800, 2000, 10000];
A_extra = [.07, .359, .935, 2.87, 53.478];
% The function Spcoef must be included from the 'all files' folder provided
% by Prof. Erik Van Vleck
[a,b,c] = Spcoef(f, A);
[d,e,f] = Spcoef(f_log, A_log);
interval=[];
fprintf('Differences between A(f) and S(f), S_log(f) respectively are:\n');
for i=1:length(A_extra)
    % The function Svalue must be included from the 'allfiles' folder provided
    % by Prof. Erik Van Vleck
    diff = A_extra(i)-Svalue(f, A, a, b, c, f_extra(i), interval);
    diff_log =log10(A_extra(i))-Svalue(f_log, A_log, d, e, f, log10(f_extra(i)), interval);
    fprintf('A(f)-S(f) is  %f at point %i\n',diff,i);
    fprintf('A(f)-S_log(f) is %f at point %i\n', diff_log, i);
end;
