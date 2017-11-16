clc
clear
% Solving by Point Gauss Siedel SOR
% To find steady state temperature distribution in a 2D slab using various
% methods
%Defining the parameters of the 2D slab

rho= 8960; %kg/m^3
Cp= 382; %J/(kg-K)
k= 400; %W/(m-K)
A=k/(rho*Cp); %m^2/s
h= 1; %m
l= 2; %m

dx=input('Enter the mesh size in the X direction : ');
dy=input('Enter the mesh size in the Y direction : ');
N=1+((l-0)/dx);
M=1+((h-0)/dy);

%Defining the temperature matrix 
Tnew=zeros(N,M);
Told=zeros(N,M);

%Defining the Boundary Conditions
Tnew=BC(Tnew);
Told=BC(Told);

%Defining the initial conditions
fprintf('\nInitial Conditions')
fprintf('\nEnter A for T(i,j)=0 degree Celcius')
fprintf('\nEnter B for T(i,j)=200 degree Celcius')
I=input('\nAnswer: ','s');
Tnew=IC(Tnew,I,dx,dy);
Told=IC(Told,I,dx,dy);

s1=(1/(dx^2));
s2=(1/(dy^2));
%Solving by Point Gauss Siedel
p=0;
a=((cos(pi/N-1)+(dx/dy)^2*cos(pi/M-1))/(1+(dx/dy)^2))^2;
beta=(2-2*((1-a)^(0.5)))/a;
beta=beta*1.3;
%relaxation variable
count=0;
while p~=1
    p=1;
    count=count+1;
    Told=Tnew;
    for i=2:(N-1)
        for j=(2:M-1)
        Tnew(i,j)=beta*(-1*Told(i+1,j)*s1-1*Told(i,j+1)*s2-1*Tnew(i-1,j)*s1-1*Tnew(i,j-1)*s2)/(-2*s1-2*s2)+(1-beta)*Told(i,j);
        end
    end
    
    for i=2:(N-1)
        for j=(2:M-1)
            if abs(Tnew(i,j)-Told(i,j))>0.001
                p=0;
            end
        end
    end
end

fprintf('The number of iterations required are %f\n',count);
pcolor(Tnew);



