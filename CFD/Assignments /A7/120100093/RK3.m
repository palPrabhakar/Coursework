function [y,t] = RK3(y_initial,x,n,h);
tic;
y(1)=y_initial;

for i=1:n
    k1=-y(i)+(x(i)^0.1)*(x(i)+1.1);
    k2=(-y(i)-0.5*h*k1)+((x(i)+0.5*h)^0.1)*(x(i)+0.5*h+1.1);
    k3=(-y(i)+h*k1-2*h*k2)+(x(i+1)^0.1)*(x(i+1)+1.1);
    y(i+1)=y(i)+(h/6)*(k1+4*k2+k3);    
end
t=toc;