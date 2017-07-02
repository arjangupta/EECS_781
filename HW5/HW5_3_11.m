% Author: Arjan Gupta
% Description: HW#5, 3.11

clear all; close all; clc;

m = [7, 10, 13];
N = [2*m(1)+1, 2*m(2)+1, 2*m(3)+1];
error_max5 = zeros(3,1);
error_max1 = zeros(3,1);
ends = [5,1]; 

% PART A of problem
% Find maximum error for P_n(x) on interval [-5,5]
fprintf('Part A: Maximum error for P_n(x) on [-5,5]\n');
for i=1:length(N)
    % Create m equally spaced different data points
    x = linspace(-5, 5, N(i));
    % Initialize an array that will hold the values of f evaluated at x
    f = zeros(1,length(x));
    for j=1:length(x)
        % Evaluate f at point x
        f(j) = 1/(1+x(j)^2);
    end;
    p=lagrangepoly(x, f);
    x1=linspace(-5,5,1001);
    for j=1:length(x1)
        diff = abs((1/(1+x1(j)^2))-polyval(p,x1(j)));
        if(diff > error_max5(i))
            error_max5(i) = diff;
        end;      
    end;
    fprintf('N = %i, error: %f\n',N(i),error_max5(i));
end;

fprintf('\n');

% PART B of problem
% Find maximum error for P_n(x) on interval [-1,1]
fprintf('Part B: Maximum error for P_n(x) on [-1,1]\n');
for i=1:length(N)
    % Create m equally spaced different data points
    x = linspace(-5, 5, N(i));
    % Initialize an array that will hold the values of f evaluated at x
    f = zeros(1,length(x));
    for j=1:length(x)
        % Evaluate f at point x
        f(j) = 1/(1+x(j)^2);
    end;
    p=lagrangepoly(x, f);
    x1=linspace(-1,1,1001);
    for j=1:length(x1)
        diff = abs((1/(1+x1(j)^2))-polyval(p,x1(j)));
        if(diff > error_max1(i))
            error_max1(i) = diff;
        end;      
    end;
    fprintf('N = %i, error: %f\n',N(i),error_max1(i));
end;