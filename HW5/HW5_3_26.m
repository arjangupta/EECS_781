% Author: Arjan Gupta
% Description: HW#5, 3.26

clear; close all; clc;

% Initialize data-set size 
n = 15;
% Assign interpolation points and values at those points
I = 605:1065;
T = [605, 645, 725, 795, 825, 845, 855, 875, 895, 905, 915, 925, 955, 1015, 1065];
C = [.622, .639, .668, .694, .73, .812, .907, 1.336, 2.169, 2.075, 1.598, 1.211, .672, .603, .601];
% Initialize extra data points
T_extra = [685, 765, 865, 885, 935, 975];
C_extra = [.655, .679, 1.044, 1.881, .916, .615];
[a,b,c] = Spcoef(T, C);
interval=[];


fprintf('Difference between C(T) and S(T) is :\n');
S1 = zeros(length(T_extra),1);
for i=1:length(T_extra)
    S1(i) = Svalue(T, C, a, b, c, T_extra(i), interval);
    diff = C_extra(i)-Svalue(T, C, a, b, c, T_extra(i), interval);
    fprintf('  %f, at point %i\n',diff,i);
end;


S2 = zeros(length(I), 1);
for j=1:length(I)
    S2(j) = Svalue(T, C, a, b, c, I(j), interval);
end;

figure
plot(I,S2)
hold on
plot(T_extra,S1,'*')