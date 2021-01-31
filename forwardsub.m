function [c] = forwardsub(L,rhs)

n = length(L);
c = zeros(n,1);
c(1) = rhs(1);

for i = 2:n
    c(i) = rhs(i)-L(i,i-1)*c(i-1);
end
end