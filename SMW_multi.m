function Ainv=SMW_multi(U,V)

[m,n]=size(U);
M=eye(n)-V*U;

Ainv=eye(m)+U*inv(M)*V;





