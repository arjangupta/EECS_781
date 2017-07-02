function [A,flag,pivots,Cond] = Factor(A)
%  Factor decomposes the matrix A using Gaussian elimination and
%  optionally estimates its condition number.  Factor is used in
%  conjunction with Solve to solve A*x = b.
%
%  Input arguments:
%     A       = matrix A of neq rows and columns to be factored.
%
%  Output arguments:
%     A       = contains the upper triangular matrix U in its upper
%               portion (by rows) and a permuted version of a lower
%               triangular matrix (I - L).  The factorization is
%               such that (permutation matrix)*A = L*U.
%    flag     = reports success or the reason for failure.  flag = 0
%               indicates success. If flag > 0, a zero pivot occurred
%               at equation number flag and the computation was
%               terminated.
%    pivots   = record of row interchanges.  The entry pivots(neq) =
%               (-1)^(number of row interchanges).
%
%  When flag > 0, the determinant of A is 0 and when flag = 0,
%               det(A) = pivots(neq) * A(1,1) *  ... * A(neq,neq).
%
%  Optional output argument:
%    Cond     = when flag >= 0, an estimate of the condition number
%               of A in the infinity norm.

[neq,cols] = size(A);
flag = 0;
pivots = zeros(neq,1);
pivots(neq) = 1;

if nargout == 4
  % Initialize Cond for A that is numerically singular.
  Cond = realmax;

  % Compute the infinity norm of A before the matrix is
  % overwritten by its factorization.
  Anorm = norm(A,inf);
end

if neq == 1 
  if A(1,1) == 0
    flag = 1;
  elseif nargout == 4
    Cond = 1;
  end
  return
end

%  Gaussian elimination with partial pivoting.
for k = 1:neq-1

  % Determine the row m containing the largest element in
  % magnitude to be used as a pivot and its magnitude biggest.
  [biggest,occurred] = max(abs(A(k:neq,k)));
  m = occurred + k - 1;

  % If all possible pivots are zero, A is numerically singular.
  if biggest == 0
    flag = k;
    return 
  end
  pivots(k) = m;
  if m ~= k
    % Interchange the current row k with the pivot row m.
    A([m k],k:neq) = A([k m],k:neq);
    pivots(neq) = - pivots(neq);
  end

  % Eliminate subdiagonal entries of column k.
  for i = k+1:neq
    t = A(i,k) / A(k,k);
    A(i,k) = - t;
    if t ~= 0
      A(i,k+1:neq) = A(i,k+1:neq) - t * A(k,k+1:neq);
    end
  end
end

if A(neq,neq) == 0
  flag = neq;
  return
end

if nargout == 4
  % Estimate the condition number of A by computing the infinity
  % norm of A directly and a lower bound for the norm of A^(-1).
  % A lower bound for the norm of A^(-1) is provided by the ratio
  % norm(y)/norm(d) for any vectors such that A*y = d and d ~= 0.
  % A "large" ratio is obtained by computing y as one iteration of
  % inverse iteration for the smallest singular value of A, i.e.,
  % by solving for y such that (A'*A)*y = e.  This exploits the
  % fact that an LU decomposition of A can be used to solve the
  % linear system A'*d = e as well as A*y = d.  The entries of e
  % are +1 or -1 with the sign chosen during the computation of d
  % to increase the size of the entry of d and so make a "large"
  % lower bound for the norm of A^(-1) more likely.

    % Solve A'*d = e using the decomposition of A.
    d = zeros(neq,1);
    d(1) = -1 / A(1,1);
    for k = 2:neq
      t = A(1:k-1,k)' * d(1:k-1);
      if t < 0
        ek = -1;
      else
        ek =  1;
      end
      d(k) = -(ek + t) / A(k,k);
    end
    for k = neq-1:-1:1
      d(k) = d(k) + A(k+1:neq,k)'*d(k+1:neq);
      m = pivots(k);
      d([m k]) = d([k m]);
    end

    % Solve A*y = d.
    y = Solve(A,pivots,d);

    % Compute the infinity norms of the vectors.
    dnorm = norm(d,inf);
    ynorm = norm(y,inf);

    Cond = max(Anorm * ynorm / dnorm, 1);

  end

