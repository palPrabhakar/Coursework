function [y,t]= for_eu(dt)
T=5;
y(1)=1; 
t(1)=0;
z(1)=0;
n=T/dt;

x=zeros(2,1);
A=zeros(2);
B=zeros(2,1);

A=[1/dt -1; 16 1/dt]; 
B=[y(1); z(1)].*(1/dt);

At=[8*dt 1; 1 -0.5*dt];
Bt=[z(1)-8*dt*y(1); y(1)+0.5*dt*z(1)];

for i=1:n
    t(i+1)=t(i)+dt;
    
    %y(i+1)=(1-.5*dt)*y(i); %Explicit
    %y(i+1)=y(i)/(1+.5*dt); % implicit
    %y(i+1)=((1-.25*dt)*y(i))/(1+.25*dt); %Trapezoidal 
   
    %{
    %RK
    k1=-0.5*y(i);
    k2=-0.5*(y(i)+(dt/2)*k1);
    k3=-0.5*(y(i)+(dt/2)*k2);
    k4=-0.5*(y(i)+(dt/2)*k3);
    y(i+1)=y(i)+((dt/6)*(k1+2*k2+2*k3+k4));
    %}
    
    %*****2nd Question******%
    
    y(i+1)=y(i)+z(i)*dt; 
    z(i+1)=z(i)-16*dt*y(i);
    t(i+1)=t(i)+dt;
    
    
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
end