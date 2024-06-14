function x=Tridiagonal_Toeplitz_Fast_Solver(alpha,beta,gamma,n,b,t)

T=toeplitz([alpha beta zeros(1,n-2)],[alpha gamma zeros(1,n-2)]);
if t==1
    T=T';
    alpha=T(1,1);
    beta=T(2,1);
    gamma=T(1,2);
end
J=[zeros(n-1,1) eye(n-1);eye(1) zeros(1,n-1)];
A_bar=J*T;
A11=A_bar(1:n-1,1:n-1);
w=A_bar(n,1:n-1);
p=A_bar(1:n-1,n);
b1=b(1);
b2=b(2:end);
v=Backward_Substitution_System_Solver_For_Tridiagonal_Toeplitz(beta,alpha,gamma,b2);
u=Backward_Substitution_System_Solver_For_Tridiagonal_Toeplitz(beta,alpha,gamma,p);
x(n)=(w*v-b1)/(w*u);
x(1:n-1)=v-x(n)*u;
x=x';
end