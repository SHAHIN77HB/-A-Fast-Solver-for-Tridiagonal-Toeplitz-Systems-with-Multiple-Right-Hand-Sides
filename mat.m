function [ A ] = mat ( x , k )
%Shahin Hasan Beigi 
%Tarbiat Modares University

%The mat(.) function takes the vecor x and return the ( N/k * k ) matrix A
%The input k is number of columns that you want.
N=length(x);
r=N/k;
for i=1:k
    A(:,i)=x((i-1)*r+1:i*r,1);
end


end

