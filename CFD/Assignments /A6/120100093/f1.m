function [y] = f1(a,b,x_max,h)

x(1)=a;
y(1)=0;

n=(x_max-x(1))/h;    

for i=1:n
    x(i+1)=x(i)+h;
    y(i+1)=x(i+1)*sin(x(i+1));
end



 