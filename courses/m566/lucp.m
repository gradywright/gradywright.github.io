function [A,p,q,gf] = lucp(A)
%lucp  Computes the LU decomposition of A with complete pivoting
%
%   [A,p,q] = lucp(A) Computes the LU decomposition of A with complete
%   pivoting using vectorization. The factors L and U are returned in the
%   output A, and the permutation of the rows and columns from complete
%   pivoting are recorded in the vectors p and q, respectively. Here L is
%   assumed to be unit lower triangular, meaning it has ones on its
%   diagonal. The resulting decomposition can be extracted from A and p as
%   follows:
%       L = eye(length(LU))+tril(LU,-1);     % L with ones on diagonal
%       U = triu(LU);                        % U 
%       P = p*ones(1,n) == ones(n,1)*(1:n);  % Permutation matrix
%   A is then given as L*U = P*A, or P'*L*U = A.
%
%   Use this function in conjuction with backsub and forsub to solve a
%   linear system Ax = b.

n = size(A,1);
p = (1:n)';
q = p;
gf = zeros(n-1,1);
maxA = max(abs(A(:)));

for k=1:n-1
    % Find the row in column k that contains the largest entry in magnitude
%     k
%     U = A;
%     for jj = 1:k-1
%         for ii = jj+1:n
%             U(ii,jj) = 0;
%         end
%     end
%     latex_rat_table(U)
    [temp,rpos] = max(abs(A(k:n,k:n)));
    [gf(k),cpos] = max(abs(temp));
    rpos = rpos(cpos);
    row2swap = k-1+rpos;
    col2swap = k-1+cpos;
    % Swap the rows of A and p (perform partial pivoting)
    A([row2swap, k],:) = A([k, row2swap],:);
    p([row2swap, k]) = p([k, row2swap]);
    A(:,[col2swap, k]) = A(:,[k, col2swap]);
    q([col2swap, k]) = q([k, col2swap]);
%     U = A;
%     for jj = 1:k-1
%         for ii = jj+1:n
%             U(ii,jj) = 0;
%         end
%     end
%     latex_rat_table(U)    
%     latex_rat_table(p)    
%     latex_rat_table(q)    
    % Perform the kth step of Gaussian elimination
    J = k+1:n;
    A(J,k) = A(J,k)/A(k,k);
    A(J,J) = A(J,J) - A(J,k)*A(k,J);
%     L = eye(size(A));
%     U = A;
%     for jj = 1:k
%         for ii = jj+1:n
%             L(ii,jj) = A(ii,jj);
%             U(ii,jj) = 0;
%         end
%     end
%     latex_rat_table(U)
%     latex_rat_table(L)    
end
gf = max(gf(:))/maxA;

gf = max(abs(A(:)))/maxA;

end
