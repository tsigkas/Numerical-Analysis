%Does it work
N=6
 %create arbitrary tridiagonal matrices
 main= N*ones(N,1) + rand(N,1);
 upper= rand(N-1,1);
 lower= rand(N-1,1);
 A = diag(main,0) + diag(upper,1) - diag(lower,-1);

%create a rhs which has the solution of the 1-vector
%compare to the solution from Thomas
known_x= ones(N,1);
manuf_b= A*known_x;
    
[L,U]=thomas(A);
c=forwardsub(L,manuf_b);
xthom=backsub(U,c)

%Cost analysis
N=500:10:1000;
ts=0*N;

for j=1:length(N)
   
    main= N(j)*ones(N(j),1) + rand(N(j),1);
    upper= rand(N(j)-1,1);
    lower= rand(N(j)-1,1);
    
    A = diag(main,0) + diag(upper,1) - diag(lower,-1);
    
    known_x= ones(N(j),1);
    manuf_b= A*known_x;
    
    %Solve 50 times and take average
    tic
    for i = 1:50
        
    [L,U]=thomas(A);
    c=forwardsub(L,manuf_b);
    xthom=backsub(U,c);
    
    end
   t=toc;
   ts(j)=t/50;
end

%regression line
linfit = polyfit(N,ts,1);
plot(N,ts,N, linfit(1)*N + linfit(2))
xlabel("Size of the tridiagonal matrix (n)")
ylabel("Time (seconds)")
legend("Time to solve the tridiagonal system", "Linear regression line", "location", "northwest")













