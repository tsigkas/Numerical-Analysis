p = @(x) -4./x;
q = @(x) -2./x.^2;
r = @(x) 2*log(x)./(x.^2);
a = 1; %Lower & upper bound
b = 2;
y1 = 1/2; %Boundary values
y2 = log(2);
n = linspace(3,40,38); %Values of n that we solve for

y_ex_f = @(x) (4./x - 2./(x.^2) + log(x) -3/2); %Exact solution

e = zeros(length(n),1); %error, inf norm
er = zeros(length(n),1); %error ratio
for i = 1:length(n)
    [A, xVals, rhs] = MVrep(p, q, r, y1, y2, a, b, n(i));
    y_exVals = y_ex_f(xVals);
    y_aprVals = zeros(n(i)+1,1);
    y_aprVals(1) = y1;
    y_aprVals(end) = y2;
    %LU-factorize A with Thomas
    [L, U] = thomas(A);
    %Forward, backward (Lc=rhs, UyVals=c)
    c = forwardsub(L, rhs);
    y_aprVals(2:n(i)) = backsub(U,c);
    e(i) = norm(y_exVals - y_aprVals, Inf);
    plot(xVals, y_aprVals)
    hold on
    if i > 1
        er(i) = e(i-1)/e(i);
    end
end

x100 = linspace(a,b);
plot(x100, y_ex_f(x100), "Linewidth", 1.5) %Add the real function to the plot
xlabel("x");
ylabel("y");
legend("n= 3","n= 5","n= 10", "exact solution", "location", "southeast")

plot(1./n.^2,e)
xlabel("1/n^2");
ylabel("Inf norm");
