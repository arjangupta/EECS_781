function x = Solve(A,pivots,b)
%  Solve solves A*x = b, a system of neq linear equations in neq
%  unknowns using the decomposition obtained from a successful call
%  to Factor.
%
%  Input arguments:
%     A           = output of Factor.  Triangular decomposition
%                   of the coefficient matrix.
%     pivots      = output of Factor. Record of row interchanges.
%     b           = right hand side vector b of length neq.
%
%  Output argument:
%     x           = solution vector of the same size as b.

neq = length(b);
x = b;

if neq == 1

  x(1) = x(1)/A(1,1);

else

  %  Forward elimination.
  for k = 1:neq-1
    m = pivots(k);
    x([m k]) = x([k m]);
    x(k+1:neq) = x(k+1:neq) + A(k+1:neq,k)*x(k);
  end

  %  Back substitution.
  x(neq) = x(neq) / A(neq,neq);
  for i = neq-1:-1:1
    x(i) = (x(i) - A(i,i+1:neq)*x(i+1:neq)) / A(i,i);
  end

end

