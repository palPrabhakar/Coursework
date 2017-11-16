function [y,t] = RK2(y_initial,x,n,h);

tic;
y(1)=y_initial;

for i=1:n
    k1=-y(i)+(x(i)^0.1)*(x(i)+1.1);
    k2=(-y(i)-h*k1)+(x(i+1)^0.1)*(x(i+1)+1.1);
    y(i+1)=y(i)+h*(0.5*k1+0.5*k2);    
end

t=toc;

