function [Q,H] = arnoldi3(A,b,n)
% ARNOLDI   Arnoldi iteration for Krylov subspaces.
% Input:
%   A    square matrix (m by m)
%   b    initial vector
%   n    number of iterations
% Output: 
%   Q    orthonormal basis of Krylov space (m by n+1)
%   H    upper Hessenberg matrix, A*Q(:,1:n)=Q*H (n+1 by n)

m = length(A);
Q = zeros(m,n+1);  H = zeros(n+1,n)
Q(:,1) = b/norm(b);
for j = 1:n
  v = A*Q(:,j);
  for i = 1:j
    H(i,j) = Q(:,i)'*v;
    v = v - H(i,j)*Q(:,i);
  end
  H(j+1,j) = norm(v);
  Q(:,j+1) = v/H(j+1,j);
end

