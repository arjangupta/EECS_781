function [b,c,d] = Spcoef(x,f)
%  Calculate coefficients defining a smooth cubic interpolatory spline S.
%
%  Input arguments:
%     x   = vector of n values of the independent variable ordered
%           so that  x(1) < x(2) < ... < x(n).
%     f   = vector of values of the dependent variable.
%
%  Output arguments:
%     b   = column vector of S'(x(i)) values.
%     c   = column vector of S''(x(i))/2 values.
%     d   = column vector of S'''(x(i)+)/6 values (i < n).

x = x(:);
f = f(:);
n = length(x);
if length(f) < n
  error('f has fewer entries than x.')
end
b = zeros(n,1);
c = zeros(n,1);
d = zeros(n,1);

h = diff(x);
if any(h <= 0)
  error('Spcoef requires that x(1) < x(2) < ... < x(n).')
end

%  Calculate coefficients for the tridiagonal system: store
%  sub-diagonal in b, diagonal in d, difference quotient in c.

b(1:n-1) = h;
diff_f = diff(f);
c(1:n-1) = diff_f ./ h;
if n == 2
  b(1) = c(1);
  c(1) = 0;
  d(1) = 0;
  b(2) = b(1);
  c(2) = 0;
  return
end
d(1) = 2 * b(1);
for i = 2:n-1
  d(i) = 2 * (b(i) + b(i-1));
end
d(n) = 2 * b(n-1);

%  Calculate estimates for the end slopes using polynomials
%  that interpolate the data nearest the end.
fp1 = c(1) - b(1) * (c(2) - c(1)) / (b(1) + b(2));
if n > 3
  fp1 = fp1 + b(1) * ((b(1) + b(2)) * (c(3) - c(2)) / ...
       (b(2) + b(3)) - c(2) + c(1)) / (x(4) - x(1));
end      
          
fpn = c(n-1) + b(n-1) * (c(n-1) - c(n-2)) / (b(n-2) + b(n-1));
if n > 3
  fpn = fpn + b(n-1) * (c(n-1) - c(n-2) - (b(n-2) + b(n-1))*...
       (c(n-2) - c(n-3)) / (b(n-2) + b(n-3))) / (x(n) - x(n-3));
end

%  Calculate the right hand side and store it in c.
c(n) = 3 * (fpn - c(n-1));
for i = n-1:-1:2
  c(i) = 3 * (c(i) - c(i-1));
end
c(1) = 3 * (c(1) - fp1);

%  Solve the tridiagonal system.
for k = 2:n
  p = b(k-1) / d(k-1);
  d(k) = d(k) - p * b(k-1);
  c(k) = c(k) - p * c(k-1);
end
c(n) = c(n) / d(n);
for k = n-1:-1:1
  c(k) = (c(k) - b(k) * c(k+1)) / d(k);
end

%  Calculate the coefficients defining the spline.
d(1:n-1) = diff(c) ./ (3 * h);
b(1:n-1) = diff_f ./ h - h .* (c(1:n-1) + h .* d(1:n-1));
b(n) = b(n-1) + h(n-1) * (2 * c(n-1) + h(n-1) * 3 * d(n-1));

