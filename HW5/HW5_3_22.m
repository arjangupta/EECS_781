% Author: Arjan Gupta
% Description: HW#5, 3.22

clear; close all; clc;

% Declare data-set size
n = 11;
% Assign interpolation points and values at those points
T = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
P = [.006107,.012277,.023378,.042433,.073774,.12338,.19924, 31166,.47364,.70112,1.01325];
% Initialize extra data points
T_extra = [5, 45, 95];
P_extra = [.008721, .095848, .84528]; 
% The function Spcoef must be included from the 'all files' folder provided
% by Prof. Erik Van Vleck
[a,b,c] = Spcoef(T,P); 
interval=[];
fprintf('Difference between P(x) and S(x) is:\n');
for i=1:length(T_extra)
    % The function Svalue must be included from the 'allfiles' folder provided
    % by Prof. Erik Van Vleck
    diff = P_extra(i)-Svalue(T, P, a, b, c, T_extra(i), interval);
    fprintf('  %f at point %i\n',diff,i);
end;
