function [X,time]=Tri_Toeplitz_MultiRHS_Fast_Solver(T,B)
[n,m]=size(B);
A=kron(eye(size(B,2)),T);
b=vec(B);
[nA,~]=size(A);
E1=[zeros(n-1,1) zeros(n-1); T(1,2) zeros(1,n-1)];
E2=[zeros(1,n-1) T(2,1) ; zeros(n-1) zeros(n-1,1)];
J=[zeros(m-1,1) eye(m-1);0 zeros(1,m-1)];
E=kron(J,E1)+kron(J',E2);
A_hat=E+A;
%%
t1=T(1,1);t2=T(2,1);t3=T(1,2);
for i=1:m-1
    alpha(:,i)=alpha_beta_e(i*n,nA,t3);
    beta(:,i)=alpha_beta_e(i*n+1,nA,t2);
    I_alpha(:,i)=alpha_beta_e(i*n+1,nA,1);
    I_beta(:,i)=alpha_beta_e(i*n,nA,1);
end
tic
for j=1:2*m-2
    if j<m
        v(:,j)=Tridiagonal_Toeplitz_Fast_Solver(t1,t2,t3,nA,I_alpha(:,j),1);
    else
        v(:,j)=Tridiagonal_Toeplitz_Fast_Solver(t1,t2,t3,nA,I_beta(:,j-m+1),1);
    end
end
u=[alpha beta];
%%

Q_inv=SMW_multi(u,v');
x=Tridiagonal_Toeplitz_Fast_Solver(t1,t2,t3,nA,Q_inv*b,0);
X=mat(x,m);
time=toc;