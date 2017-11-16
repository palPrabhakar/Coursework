clc
clear
% To find steady state temperature distribution in a 2D slab using line Gauss
% Siedel SOR Method
%Defining the parameters of the 2D slab

rho= 8960; %kg/m^3
cp= 382; %J/(kg-K)
k= 400; %W/(m-K)
alpha=k/(rho*cp); %m^2/s
H= 1; %m
L= 2; %m
beta=1.21;

dx=input('Enter the mesh size in the X direction : ');
dy=input('Enter the mesh size in the Y direction : ');
n=1+((L-0)/dx);
m=1+((H-0)/dy);

T=zeros(n,m);
T2=zeros(n,m);
T3=zeros(n,m);

%Defining the Boundary Conditions 
T=BC(T);

%Defining the initial conditions
fprintf('\nInitial Conditions')
fprintf('\nEnter A for T(i,j)=0 degree Celcius')
fprintf('\nEnter B for T(i,j)=200 degree Celcius')
I=input('\nAnswer: ','s');
T=IC(T,I,dx,dy);

f=0;                % Flag variable  
count=0;            % Counter

while f~=1
    f=1;
    T2=T;           % Previous iteration temperatures
    for j=2:m-1
        
        % Block for defining the coefficients for the line gauss method with SOR
        % Initializing the coefficients 
        
        a1=zeros(n-2,1);
        a2=zeros(n-2,1);
        a3=zeros(n-2,1);
        b=zeros(n-2,1);
        
        % Writing the coefficients of the tridiagonal system of dfference
        % equation obtained for line gauss method with SOR
        
        for i=2:(n-3)
            a1(i,1)=-beta;
            a2(i,1)=4;
            a3(i,1)=-beta;
            b(i,1)=beta*(T(i+1,j-1)+T(i+1,j+1)) + (1-beta)*4*T(i+1,j);
        end
        
        % Coefficients for the boundary points
        a1(n-2,1)=-beta;
        a3(1,1)=-beta;
        a2(1,1)=4;
        a2(n-2,1)=4;
        b(1,1)= beta*(T(2,j-1) + T(2,j+1))+beta*T(1,j) + (1-beta)*4*T(2,j);
        b(n-2,1)= beta*(T(n-1,j-1) + T(n-1,j+1)) + beta*T(n,j) + (1-beta)*4*T(n-1,j);
        
        % Solving the tridiagonalsystem using the function developed earlier
        x=tridiagonal(n-2,a1,a2,a3,b);
        
        % Writing the solution for the latest time step
        for i=2:n-1
            T(i,j)=x(1,i-1);      % Updated Temperature
        end                  
    end
    
    for i=2:(n-1)
        for j=2:(m-1)
            if abs(T2(i,j)-T(i,j))>0.001    % Convergence criteria
                f=0;
            end
        end
    end
    count =count +1;            % Counts iterations
  
end
fprintf('The number of iterations required are %f\n',count);
pcolor(T)