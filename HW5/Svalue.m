function [S,interval,flag] = Svalue(x,f,b,c,d,t,interval)
%  Evaluate the spline S at t using coefficients from Spcoef.
%
%  Input arguments:
%     x, f, b, c, d are defined as in Spcoef.
%     t             point where spline is to be evaluated.
%     interval      On the first call to Svalue, set to [].
%                   Thereafter use the value returned from
%                   the last call.
%
%  Output arguments:
%     S        = value of spline at t.
%     interval = index satisfying  x(interval) <= t < x(interval+1)
%                unless t is outside data interval (see flag).
%     flag     =  0  normal return;
%              =  1  t < x(1).
%              =  2  t > x(n).

n = length(x);

if isempty(interval) 
  last_interval = round((n + 1) / 2);
else
  last_interval = interval;
end

%  Search for correct interval for t.
if     (t < x(1))
  flag = 1;
  interval = 1;
elseif (t < x(n))                 % x(1) <= t < x(n).
  flag = 0;
  if t >= x(last_interval)
    for j = last_interval:n-1
      if t < x(j+1)
        interval = j;
        break
      end
    end
  else
    for j = last_interval-1:-1:1
      if t >= x(j)
        interval = j;
        break
      end
    end
  end
elseif (t > x(n))
  flag = 2;
  interval = n-1;
else                              % t = x(n)
  flag = 0;
  interval = n-1;
end

%  Evaluate cubic polynomial.
dt = t - x(interval);
S = f(interval) + dt * (b(interval) + dt * (c(interval) + ...
                  dt * d(interval)));



