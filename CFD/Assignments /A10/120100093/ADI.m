
clc
clear
% To find steady state temperature distribution in a 2D slab ADI
%Defining the parameters of the 2D slab

rho= 8960; %kg/m^3
cp= 382; %J/(kg-K)
k= 400; %W/(m-K)
alpha=k/(rho*cp); %m^2/s
H= 1; %m
L= 2; %m

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
    T2=T;           % Initializing intermediate temperature T2
    for j=2:m-1
        
        % Block for defining the coefficients for the ADI method
        % Initializing the coefficients 
        
        a1=zeros(n-2,1);
        a2=zeros(n-2,1);
        a3=zeros(n-2,1);
        b=zeros(n-2,1);
        
        % Writing the coefficients of the tridiagonal system of dfference
        % equation obtained for ADI method
        
        for i=2:(n-3)
            a1(i,1)=-1;
            a2(i,1)=4;
            a3(i,1)=-1;
            b(i,1)=T2(i+1,j-1)+T(i+1,j+1);
        end
        
        % Coefficients for the boundary points
        a1(n-2,1)=-1;
        a3(1,1)=-1;
        a2(1,1)=4;
        a2(n-2,1)=4;
        b(1,1)=T2(2,j-1)+T(2,j+1)+T2(1,j);
        b(n-2,1)=T2(n-1,j-1)+T(n-1,j+1)+T2(n,j);
        
        % Solving the tridiagonalsystem using the function developed earlier
        x=tridiagonal(n-2,a1,a2,a3,b);
        
        % Writing the solution for the latest time step
        for i=2:n-1
            T2(i,j)=x(1,i-1);      % T2 is the intermediate temp at the count + 1/2 time step
        end                  
    end
    
    T3=T;               % Initializing final temperature T3
    for i=2:n-1
        
        % Block for defining the coefficients for the ADI method
        % Initializing the coefficients 
        
        a1=zeros(m-2,1);
        a2=zeros(m-2,1);
        a3=zeros(m-2,1);
        b=zeros(m-2,1);
        
        % Writing the coefficients of the tridiagonal system of dfference
        % equation obtained for ADI method
        
        for j=2:(m-3)
            a1(j,1)=-1;
            a2(j,1)=4;
            a3(j,1)=-1;
            b(j,1)=T3(i-1,j+1)+T2(i+1,j+1);
        end
        
        % Coefficients for the boundary points
        a1(m-2,1)=-1;
        a3(1,1)=-1;
        a2(1,1)=4;
        a2(m-2,1)=4;
        b(1,1)=T3(i-1,2)+T2(i+1,2)+T3(i,1);
        b(m-2,1)=T3(i-1,m-1)+T2(i+1,m-1)+T3(i,m);
        
        % Solving the tridiagonalsystem using the function developed earlier
        x2=tridiagonal(m-2,a1,a2,a3,b);
        
        % Writing the solution for the latest time step
        for j=2:m-1
            T3(i,j)=x2(1,j-1);      % T3 is the final temp at the count + 1 time step
        end                  
    end
    
    for i=2:(n-1)
        for j=2:(m-1)
            if abs(T3(i,j)-T(i,j))>0.001    % Convergence criteria
                f=0;
            end
        end
    end
    count =count +1;            % Counts iterations
    T=T3;
   
end
fprintf('The number of iterations required are %f\n',count);
pcolor(T)