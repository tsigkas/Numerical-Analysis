function [L,U] = thomas(A)
% LU = A

     n  = length(A);
   for j = 1:n-1
      A(j+1,j) = A(j+1,j)/A(j,j); %L-component
      A(j+1,j+1) = A(j+1,j+1) - A(j+1,j)*A(j,j+1); %diagonal U-component
   end

   U = triu(A);
   L = eye(n) + tril(A,-1);
end