function [ z ] = vec ( a )

[n m]=size(a);

 z=a(:,1);
for i=2:m
   
    z=[z;a(:,i)];

end

