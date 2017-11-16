function [x]=Tri_Diagonal(A,n,b)
% Tridiagonal-Matrix Solution

x=zeros(n,1);

for i=1
    p=A(i,i);
    A(1,:)=A(1,:)/p;
    b(i,1)=b(i,1)/p;
end

for i=2:n
    q=A(i,i-1);
    A(i,:)=A(i,:)-(q*A(i-1,:));
    b(i,1)=b(i,1)-(q*b(i-1,1));
    r=A(i,i);
    A(i,:)=A(i,:)/r;
    b(i,1)=b(i,1)/r;
end

x(n,1)=b(n,1);

for i=n-1:-1:1
    x(i,1)=b(i,1)-(A(i,i+1)*x(i+1,1));
end