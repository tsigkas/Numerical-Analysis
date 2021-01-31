function [A, xVals, rhs] = MVrep(p, q, r, y0, yn, a, b, n)

% y'' = p(x)*y' + q(x)*y + r(x) (p,q,r passed as anonymous functions)
% [a,b] is the interval, n are the number of steps
% y0 and yn are the values at the boundaries

% A: (n-1)x(n-1) matrix
%xVals: (n+1)x1 vector
% b: (n-1)x1 vector

h = (b-a)/(n);

xVals = zeros(n+1,1);
xVals(1) = a;
for i = 2:n+1
    xVals(i) = xVals(i-1) + h;
end

pVals = p(xVals);
qVals = q(xVals);
rVals = r(xVals);

rhs = zeros(n-1,1);

rhs = h^2*rVals(2:n); %h^2 * r at  interior points
rhs(1) = rhs(1) - (1+pVals(2)*h/2)*y0; %lower bound
rhs(end) = rhs(end) - (1-pVals(n)*h/2)*yn; % upper bound

main_diag = -qVals(2:n)*h^2 -2*ones(n-1,1);
lower_diag = pVals(3:n)*h/2 + ones(n-2,1);
upper_diag = -pVals(2:n-1)*h/2 + ones(n-2,1);

A = diag(main_diag,0) + diag(upper_diag,1) + diag(lower_diag,-1);