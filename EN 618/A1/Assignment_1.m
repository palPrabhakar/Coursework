clear all
clc

t(1)=0;
y(1)=1;
dt=[1 0.1 .01 .001 .0001]; 
T=5;
n=T/dt(5);
y_exact(1)=exp(-0.5*t(1));
for i=1:n
    %y1(i+1)=(1-.5*dt1)*y1(i); %explicit

    
    %y(i+1)=y(i)/(1+.5*dt); % implicit
    %y(i+1)=((1-.25*dt)*y(i))/(1+.25*dt); %Trapezoidal 
    
    % RK 4
    %k1=-0.5*y(i);
    %k2=-0.5*(y(i)+(dt/2)*k1);
    %k3=-0.5*(y(i)+(dt/2)*k2);
    %k4=-0.5*(y(i)+(dt/2)*k3);
    %y(i+1)=y(i)+((dt/6)*(k1+2*k2+2*k3+k4));
    
    t(i+1)=t(i)+dt(5);
    y_exact(i+1)=exp(-0.5*t(i+1));
end

[y1,t1]=for_eu(dt(1));
[y2,t2]=for_eu(dt(2));
[y3,t3]=for_eu(dt(3));
[y4,t4]=for_eu(dt(4));
[y5,t5]=for_eu(dt(5));


plot(t1,y1,'r-o');
hold on
plot(t2,y2,'g');
plot(t3,y3,'y');
plot(t4,y4,'m');
plot(t5,y5,'k');
plot(t5,y_exact,'b');
legend('dt=1','dt=0.1','dt=0.01','dt=0.001','dt=0.0001','y_exact')
xlabel('time') % x-axis label
ylabel('y values') % y-axis label
title('Trapezoidal')


