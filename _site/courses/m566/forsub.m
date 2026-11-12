function x = forsub(L,b)
%forsub  Solves the unit lower triangular system Lx=b using forward
%        substitution.
%
%   x = forsub(L,b) Solves the unit lower triangular system Lx=b using
%   forward substition.
%
%   Use in conjuction with lupp and backsub to solve a general system Ax=b

n = size(L,1);   % How big is the system
x = zeros(n,1);  % Initial the solution vector

x(1) = b(1);
for j = 2:n
    x(j) = b(j) - L(j,1:j-1)*x(1:j-1);
end
end


