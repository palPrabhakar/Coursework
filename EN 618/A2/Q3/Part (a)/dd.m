clc
clear all
hold on
%Explicit Euler
alpha=0.0001;L=5;x0=0.5; sigma=0.01;u=0.2;
n=2001;
dx=L/(n-1)
x=0:dx:L;
T=exp(-((x-x0).^2)/sigma);
Tn=T;
Tn(1)=0; Tn(n)=0;
Tnp=Tn;
%dt=(dx.^2)/(4*alpha)
dt=0.002
N=floor(8/dt)+1
while j<N
    for i=2:n-1
        Tnp(i)=Tn(i)+dt*(alpha*[(Tn(i+1)-2*Tn(i)+Tn(i-1))/(dx.^2)]-[u*(Tn(i)-Tn(i-1))/(dx)]);
    end 
    Tn=Tnp;
    j=j+1;
    Tnp;
end 
plot (x,Tnp,'green--')

clear all
%Implicit Euler
alpha=0.0001;L=5;x0=0.5; sigma=0.01;u=0.2;
n=2001;
dx=L/(n-1);
x=0:dx:L;
T=exp(-((x-x0).^2)/sigma);
Tn=T;
Tn(1)=0; Tn(n)=0;
Tnp=Tn;
Tnp(1)=0; Tnp(n)=0;
dt=0.002
N=(8/dt)+1
while j<N
    err=0;
    while (err<10^(-6))
    for i=2:n-1
        Tnp(i)=[Tn(i)/dt+(alpha/(dx.^2))*Tn(i+1)+(alpha/(dx.^2)+u/(dx))*Tnp(i-1)]/[(1/dt)+(2*alpha/(dx.^2))+u/dx];
    end
    for i=1:n
    err=err+(Tnp(i)-Tn(i)).^2;
    end 
    err=sqrt(err);
    end 
    Tn=Tnp;
    j=j+1;
    Tnp;
end 
plot (x,Tnp,'black--')

clear all
%Trapezoidal
alpha=0.0001;L=5;x0=0.5; sigma=0.01;u=0.2;
n=2001;
dx=L/(n-1);
x=0:dx:L;
T=exp(-((x-x0).^2)/sigma);
Tn=T;
Tn(1)=0; Tn(n)=0;
Tnp=Tn;
Tnp(1)=0; Tnp(n)=0;
dt=0.002
N=(8/dt)+1
while j<N
    err=0;
    while (err<10^(-6))
    for i=2:n-1
        Tnp(i)=[(1/dt-alpha/(dx.^2)-u/(2*dx))*Tn(i)+(0.5*alpha/(dx.^2))*(2*Tn(i+1))+(0.5*alpha/(dx.^2)+u/(2*dx))*(Tnp(i-1)+Tn(i-1))]/[(1/dt)+(alpha/(dx.^2))+u/(2*dx)];
    end
    for i=1:n
    err=err+(Tnp(i)-Tn(i)).^2;
    end 
    err=sqrt(err);
    end 
    Tn=Tnp;
    j=j+1;
    Tnp;
end 
plot (x,Tnp,'red--')

clear all
%RK4 
alpha=0.0001;L=5;x0=0.5; sigma=0.01;u=0.2;
n=2001;
dx=L/(n-1)
x=0:dx:L;
T=exp(-((x-x0).^2)/sigma);
Tn=T;
Tn(1)=0; Tn(n)=0;
Tnp=Tn;
Tnp(1)=0; Tnp(n)=0;
dt_prime=2.79*(dx.^2)/(4*alpha);
%dt=dt_prime/2
dt=0.002
N=floor(8/dt)+1
k1(1:10000)=0;
k2(1:10000)=0;
k3(1:10000)=0;
k4(1:10000)=0;
while j<N
    for i=2:n-1
        k1(i)=dt*(alpha*(Tn(i+1)-2*Tn(i)+Tn(i-1))/(dx.^2)-u*(Tn(i)-Tn(i-1))/(dx));
    end
    for i=2:n-1
        k2(i)=dt*(alpha*([Tn(i+1)+k1(i+1)/2]-2*[Tn(i)+k1(i)/2]+[Tn(i-1)+k1(i-1)/2])/(dx.^2)-u*([Tn(i)+k1(i)/2]-[Tn(i-1)+k1(i-1)/2])/(dx));
    end 
    for i=2:n-1
        k3(i)=dt*(alpha*([Tn(i+1)+k2(i+1)/2]-2*[Tn(i)+k2(i)/2]+[Tn(i-1)+k2(i-1)/2])/(dx.^2)-u*([Tn(i)+k2(i)/2]-[Tn(i-1)+k2(i-1)/2])/(dx)); 
    end
    for i=2:n-1
        k4(i)=dt*(alpha*([Tn(i+1)+k3(i+1)/2]-2*[Tn(i)+k3(i)/2]+[Tn(i-1)+k3(i-1)/2])/(dx.^2)-u*([Tn(i)+k3(i)/2]-[Tn(i-1)+k3(i-1)/2])/(dx)); 
    end 
    for i=2:n-1
        Tnp(i)=Tn(i)+k1(i)/6+(k2(i)+k3(i))/3+k4(i)/6;
    end 
    Tn=Tnp;
    j=j+1;
    Tnp;
end 
plot (x,Tnp,'yellow--')

T2=(1/sqrt(1+4*alpha*8/sigma))*exp(-((x-0.5-u*8).^2)/(sigma+4*alpha*8));
hold on
plot (x,T2,'blue')

title ('Q3b for dx=0.025 and dt=0.025')
legend('Explicit','Implicit','Trapezoidal','RK4','Exact')
%legend('Implicit','Trapezoidal','RK4','Exact')
