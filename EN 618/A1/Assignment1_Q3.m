clear all
clc

t(1)=0;
yl(1)=0.5;
ynl(1)=0.5;
xl(1)=0.5;
xnl(1)=0.5;
dt=[.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00005 ];
T=1;

for j=1:7
    n(j)=T/dt(j);
    
    % Immplicit Euler
    for i=1:n(j)

        %Linearized
        yl(i+1)=yl(i)+(dt(j)*(yl(i)*(yl(i)-1)+dt(j)*(2*yl(i)-1)));

        %Non linearized
        ynl(i+1)=((dt(j)+1)-sqrt(dt(j)*dt(j)+2*dt(j)-4*dt(j)*ynl(i)+1))/(2*dt(j));

    end

    % Trapezoidal Rule 
    for i=1:n(j)

        %Linearized 
        xl(i+1)=xl(i)+(dt(j)*xl(i)*(xl(i)-1))/(1-dt(j)*(xl(i)-0.5));

        %Non Linearized 
        xnl(i+1)=(((2/dt(j))+1)-sqrt(((2/dt(j))+1)^2-(4*((2/dt(j))*xnl(i)+xnl(i)*(xnl(i)-1)))))/2;

    end

    y_exact=(1/(exp(1)+1));
    
    err1(j)=abs(y_exact-yl(n(j)+1));
    err2(j)=abs(y_exact-ynl(n(j)+1));
    err3(j)=abs(y_exact-xl(n(j)+1));
    err4(j)=abs(y_exact-xnl(n(j)+1));
    
end

hold on
plot(log(n),log(err1),'r');
plot(log(n),log(err2),'b');
plot(log(n),log(err3),'g');
plot(log(n),log(err4),'m');
legend('Lin. Implicit Euler','Non Lin. Implicit Euler','Lin. Trapezoidal','Non Lin. Trapezoidal')
xlabel('log(n)') % x-axis label
ylabel('log(error)') % y-axis label
title('Q3 Example 4.4')


