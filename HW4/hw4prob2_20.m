% Author: Arjan Gupta
% Description: Code for HW#4, Problem 2.20 a), MATH 781

clc, clear, close all;

a = 482317;
b = 2196.05;
y = 6708.43;

A1 = [a 0 0 0  b -b
      0 a 0 -b  0 -b
      0 0 a  b  b 0
      0 -b b y  0 0
      b 0  b 0  y 0
      -b -b 0 0 0 y];
  
%  Factor A1 and check its condition number.

[A1,flag,pivot_index,Cond] = Factor(A1);
if flag > 0
  fprintf('A1 has a zero pivot at %i\n',flag);
else
  fprintf('Cond = %e\n',Cond);

% Define the first right hand side vector B and solve the system.
  b = [15 0 -15 0 25 0]';
  x = Solve(A1,pivot_index,b);
  disp('Solution of the first system')
  disp(x')
end