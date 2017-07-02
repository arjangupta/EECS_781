% Author: Arjan Gupta
% Description: Code for HW#4, Problem 3.4 b), MATH 781

clc, clear, close all;

n = 6;

% interpolating points:
I1 = [10,20,40,60,80,100];
I2 = [5,15,30,50,70,90];

%  Assign initial values to the entries of A1 and A2.
A1 = zeros(n);
A2 = zeros(n);
 
% Fill the rest of A1 ans A2:

for i = 1:n
    for j = 1:n
        A1(i,j) = [I1(1,i)]^(j-1); 
    end
end
 
%disp(A)

for i = 1:n
    for j = 1:n
        A2(i,j) = [I2(1,i)]^(j-1); 
    end
end

%  Factor A1 and check its condition number.

[A1,flag,pivot_index,Cond] = Factor(A1);
if flag > 0
  fprintf('A1 has a zero pivot at %i\n',flag);
else
  fprintf('Cond = %e\n',Cond);

% Define the first right hand side vector B and solve the system.
  b = [1.498,2.138,2.840,2.542, 1.877, 1.201]';
  x = Solve(A1,pivot_index,b);
  disp('Solution of the first system')
  disp(x')
end

%  Factor A2 and check its condition number.

[A2,flag,pivot_index,Cond] = Factor(A2);
if flag > 0
  fprintf('A2 has a zero pivot at %i\n',flag);
else
  fprintf('Cond = %e\n',Cond);

% Define the first right hand side vector B and solve the system.
  b = [1.226,1.822,2.662,2.807, 2.210, 1.539]';
  x = Solve(A2,pivot_index,b);
  disp('Solution of the second system')
  disp(x')
end
