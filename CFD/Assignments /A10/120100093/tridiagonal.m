function [x]=tridiagonal(n,a1,a2,a3,b);
n1=size(a1);
n2=size(a2);
n3=size(a3);
if n1(1)==n
    if n2(1)==n
        if n3(1)==n
a=zeros(n,n);
a(1,1)=a2(1);a(1,2)=a3(1);a(n,n)=a2(n);a(n,n-1)=a1(n);
for i=2:1:n-1
    a(i,i-1)=a1(i);
    a(i,i)=a2(i);
    a(i,i+1)=a3(i);
end
for i=2:1:n
    a(i,i)=a(i,i)-((a(i,i-1)/a(i-1,i-1))*a(i-1,i));
    b(i)=b(i)-((a(i,i-1)/a(i-1,i-1))*b(i-1));
end
x(n)=b(n)/a(n,n);
for i=n-1:-1:1
    x(i)=(b(i)-(a(i,i+1)*x(i+1)))/a(i,i);
end
        else
        end
    else    
    end
else    
end