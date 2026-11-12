function x = backsub(U,b)
%forsub  Solves the upper triangular system Ux=b using backward
%        substitution.
%
%   x = backsub(U,b) Solves the upper triangular system Ux=b using
%   backward substition.
%
%   Use in conjuction with lupp and forsub to solve a general system Ax=b

n = size(U,1);   % How big is the system
x = zeros(n,1);  % Initial the solution vector

x(n) = b(n)/U(n,n);
for j = n-1:-1:1
    x(j) = ( b(j) - U(j,j+1:n)*x(j+1:n) )/U(j,j);
end

end
