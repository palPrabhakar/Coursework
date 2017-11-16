function [y,t] = AB4(y_initial,x,n,h)

tic;
y(1) = y_initial;;


for i = 1 : 3 % Inital values using Runge-Kutta of order 4
    k1=-y(i)+(x(i)^0.1)*(x(i)+1.1);
    k2=(-y(i)-0.5*h*k1)+((x(i)+0.5*h)^0.1)*(x(i)+0.5*h+1.1);
    k3=(-y(i)-0.5*h*k2)+((x(i)+0.5*h)^0.1)*(x(i)+0.5*h+1.1);
    k4=(-y(i)-h*k3)+(x(i+1)^0.1)*(x(i+1)+1.1);
    y(i+1)=y(i)+(h/6)*(k1+2*k2+2*k3+k4); 
end;

for i=4:n % main phase
    y(i+1) = y(i) + (h/24)*(55 *(-y(i)+(x(i)^0.1)*(x(i)+1.1))-59*(-y(i-1)+(x(i-1)^0.1)*(x(i-1)+1.1))+37*(-y(i-2)+(x(i-2)^0.1)*(x(i-2)+1.1))-9*(-y(i-3)+(x(i-3)^0.1)*(x(i-3)+1.1)));
end;
t=toc;