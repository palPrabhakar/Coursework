clear all
clc

t(1)=0;

dt=[0.1 , .01, .001, 0.0001, .0005];
T=5;


for j=1:5
    n=T/dt(j);
    y=zeros(5,n+1);
    t=zeros(5,n+1);
    y_exact=zeros(5,n+1);
    y(:,1)=1;
    y_exact(j,1)=exp(-0.5*t(j,1));
    for i=1:n 
        y(j,i+1)=(1-.5*dt(j))*y(j,i); %explicit 
        %y(i+1)=y(i)/(1+.5*dt(j)); % implicit
        %y(i+1)=((1-.25*dt(j))*y(i))/(1+.25*dt(j)); %Trapezoidal 
    
         
        k1=-0.5*y(i); %RK
        k2=-0.5*(y(i)+(dt(j)/2)*k1);
        k3=-0.5*(y(i)+(dt(j)/2)*k2);
        k4=-0.5*(y(i)+(dt(j))*k3);
        y(i+1)=y(i)+((dt(j)/6)*(k1+2*k2+2*k3+k4));
        

        t(j,i+1)=t(j,i)+dt(j);
        y_exact(j,i+1)=exp(-0.5*t(j,i+1));
    end
    %err(j)=y_exact(n+1)-y(n+1);
    %plot(log(dt),log(err));
end

hold on
plot(t(1,:),y(1,:),'r');
%plot(t(2,:),y(2,:),'b');
%plot(t(3,:),y(3,:),'g');
%plot(t(4,:),y(4,:),'y');
%plot(t(5,:),y(5,:),'m');
%plot(t(5,:),y_exact(5,:),'k');
xlabel('time') % x-axis label
ylabel('Y'); % y-axis label
title('Euler Method');
legend




