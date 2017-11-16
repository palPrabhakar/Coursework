function [y1,y2,y3] = forward_euler(a,x,h,n,l1,l2)

y1=zeros(1,n+1);
y2=zeros(1,n+1);
y3=zeros(1,n+1);

y1(1)=a;
y2(1)=a;
y3(1)=a;

for i=1:n
    if i==1
        y1(1,i+1) = y1(1,i);
    else
        y1(1,i+1) = y1(1,i) + h*((x(i)*cos(x(i)))+(y1(1,i)/x(i)));
    end
    
    y2(1,i+1) = y2(1,i) + h*((l1*y2(1,i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i)))));
    y3(1,i+1) = y3(1,i) + h*((l2*y3(1,i))+(1/(1+(x(i)*x(i))))-(l2*atan(x(i))));
end
