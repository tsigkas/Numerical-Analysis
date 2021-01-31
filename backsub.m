function [yVals] = backsub(U,c)

n = length(U);
yVals = zeros(n,1);
yVals(end) = c(end)/U(n,n);

for i = n-1:-1:1
    yVals(i) = (c(i)-U(i,i+1)*yVals(i+1))/U(i,i);
end
end