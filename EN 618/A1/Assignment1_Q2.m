clear all
clc

t(1)=0;
y(1)=1;
z(1)=0;
dt=[0.1 0.01 .001 .005 .0001]; 
T=5;
n=T/dt(5);

x=zeros(2,1);
A=zeros(2);
B=zeros(2,1);

%A=[1/dt -1; 16 1/dt]; 
%B=[y(1); z(1)].*(1/dt);

%At=[8*dt 1; 1 -0.5*dt];
%Bt=[z(1)-8*dt*y(1); y(1)+0.5*dt*z(1)];

for i=1:n
    
    %{
    %Explicit
    y(i+1)=y(i)+z(i)*dt; 
    z(i+1)=z(i)-16*dt*y(i);
    t(i+1)=t(i)+dt;
    %}
    
    %{
    %Implicit
    x=(A^-1)*B;
    y(i+1)=x(1,1);
    z(i+1)=x(2,1);
    B=x.*(1/dt);  
    %}
    
    %{
    %Trapezoidal 
    x=(At^-1)*Bt;
    y(i+1)=x(1,1);
    z(i+1)=x(2,1);
    Bt=[z(i+1)-8*dt*y(i+1); y(i+1)+0.5*dt*z(i+1)];  
    %}
    
    %{
    %RK 4 Scheme
    k1=-16*y(i); 
    k2=-16*(y(i)+(dt/2)*k1);
    k3=-16*(y(i)+(dt/2)*k2);
    k4=-16*(y(i)+(dt)*k3);
    z(i+1)=z(i)+((dt/6)*(k1+2*k2+2*k3+k4));
    
    k_1=z(i); 
    k_2=(z(i)+(dt/2)*k1);
    k_3=(z(i)+(dt/2)*k2);
    k_4=(z(i)+(dt)*k3);
    y(i+1)=y(i)+((dt/6)*(k_1+2*k_2+2*k_3+k_4));
    %}
    
    t(i+1)=t(i)+dt(5);
    y_exact(i+1)=cos(4*t(i+1));
end

%plot(t,y,'r');
%hold on
%plot(t,y_exact,'b');
    
[y1,t1]=for_eu(dt(1));
[y2,t2]=for_eu(dt(2));
[y3,t3]=for_eu(dt(3));
[y4,t4]=for_eu(dt(4));
[y5,t5]=for_eu(dt(5));


plot(t1,y1,'r');
hold on
plot(t2,y2,'g');
plot(t3,y3,'y');
plot(t4,y4,'m');
plot(t5,y5,'k');
plot(t5,y_exact,'b');
legend('dt=1','dt=0.1','dt=0.01','dt=0.001','dt=0.0001','y_exact')
xlabel('time') % x-axis label
ylabel('y values') % y-axis label
title('Explicit Euler')
