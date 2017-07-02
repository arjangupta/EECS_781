% Author: Arjan Gupta
% Description: Code for HW#4, Problem 3.6, MATH 781

clc, clear, close all;

n = 12;

% interpolating points:
I1 = [5,10,15,20,30,40,50,60,70,80,90,100];

%  Assign initial values to the entries of A1 and A2.
A1 = zeros(n);
 
% Fill the rest of A1 ans A2:

for i = 1:n
    for j = 1:n
        A1(i,j) = [I1(1,i)]^(j-1); 
    end
end
 
%disp(A)

%  Factor A1 and check its condition number.

[A1,flag,pivot_index,Cond] = Factor(A1);
if flag > 0
  fprintf('A1 has a zero pivot at %i\n',flag);
else
  fprintf('Cond = %e\n',Cond);

% Define the first right hand side vector B and solve the system.
  b = [1.226,1.498,1.822,2.138,2.662,2.840,2.807,2.542,2.210,1.877,1.539,1.201]';
  x = Solve(A1,pivot_index,b);
  disp('Solution of the first system')
  disp(x')
end

figure 
plot(x)
ylabel('Values computed');
xlabel('Data points, where 1 corresponds to 5, 12 corr. to 100');