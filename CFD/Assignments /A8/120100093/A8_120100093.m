clear all 
clc

rho=4000;
Cp=8000;
k=100;
h=20;
T_top=500;
T_amb=300;
T_initial=300;

dt = input ('Enter dt: ');
dy = input ('Enter dy: ');

alpha = k/(rho*Cp) ;
r = alpha*dt/(dy^2) ;
s = h*dy/k ;

nt = (10/dt) ;
ny = (0.2/dy) ;

y=0:dy(1):0.2;
t=0:dt(1):10;

Tft1=zeros(1,ny+1);
Tft1n=zeros(1,ny+1);
Tft2=zeros(1,ny+1);
Tft2n=zeros(1,ny+1);

Tdu1=zeros(1,ny+1);
Tdu1n=zeros(1,ny+1);
Tdu2=zeros(1,ny+1);
Tdu2n=zeros(1,ny+1);

Tim1=zeros(1,ny+1);
Tim1n=zeros(1,ny+1);
Tim2=zeros(1,ny+1);
Tim2n=zeros(1,ny+1);

Tcn1=zeros(1,ny+1);
Tcn1n=zeros(1,ny+1);
Tcn2=zeros(1,ny+1);
Tcn2n=zeros(1,ny+1);

T(1,:)=T_initial;
T(ny+1)=T_top;
T(1)=T_amb;
Tn(ny+1)=T_top;
Tn(1)=T_amb;
T1(ny+1)=T_top;
T1(1)=(k*T1(2)+h*dy(1)*T_amb)/(k+h*dy(1));
Tn1(ny+1)=T_top;

%Boundary Condition 1
Tft1 = T_amb*ones(ny+1,1);
Tft1(ny+1) = T_top;
Tft1n = T_amb*ones(ny+1,1);
Tft1n(ny+1) = T_top;

Tdu1 = T_amb*ones(ny+1,1);
Tdu1(ny+1) = T_top;
Tdu1n = T_amb*ones(ny+1,1);
Tdu1n(ny+1) = T_top;

Tim1 = T_amb*ones(ny+1,1);
Tim1(ny+1) = T_top;
Tim1n = T_amb*ones(ny+1,1);
Tim1n(ny+1) = T_top;

Tcn1 = T_amb*ones(ny+1,1);
Tcn1(ny+1) = T_top;
Tcn1n = T_amb*ones(ny+1,1);
Tcn1n(ny+1) = T_top;
    
 %-----FTCS Scheme-----%
for i = 1:nt+1
    for j = 2:ny
        Tft1n(j) = Tft1(j) + (r*(Tft1(j+1)-(2*Tft1(j))+Tft1(j-1))) ;
    end
    Tft1 = Tft1n ;
end
        
 %-----Du-Fort Frenkel Scheme-----%
        
T_old = 300*ones(ny+1,1) ;
T_old(ny+1) = 500 ;

for i = 1:nt+1
    for j = 2:ny
        Tdu1n(j) = ( T_old(j) + (2*r*(Tdu1(j+1)-T_old(j)+Tdu1(j-1))) )/(1+(2*r)) ;
    end
    T_old = Tdu1 ;
    Tdu = Tdu1n ;
end
        
%-----Implicit Scheme-----%

A = zeros(ny-1) ;
b = zeros(ny-1,1) ;

for i = 1:nt+1

    for j = 1:ny-1
        A(j,j) = -(1+(2*r)) ;
    end

    for j = 1:ny-2
        A(j,j+1) = r ;
    end

    for j = 2:ny-1
        A(j,j-1) = r ;
    end

    b(1) = -Tim1(2)-(Tim1(1)*r) ;
    b(ny-1) = -Tim1(ny)-(Tim1(ny+1)*r) ;

    for j = 2:ny-2
        b(j) = -T(j+1) ;
    end

    [X] = Tri_Diagonal(A,ny-1,b) ;

    for j = 1:ny-1
        Tcn1n(j+1) = X(j) ;
    end

    T = Tcn1n ;

end

%-----Crank-Nicholson Scheme-----%

A = zeros(ny-1) ;
b = zeros(ny-1,1) ;

for i = 1:100

    for j = 1:ny-1
        A(j,j) = 2*(1+r) ;
    end

    for j = 1:ny-2
        A(j,j+1) = -r ;
    end

    for j = 2:ny-1
        A(j,j-1) = -r ;
    end

    b(1) = (r*Tcn1(3))+(2*(1-r)*Tcn1(2)+(2*r*Tcn1(1))) ;
    b(ny-1) = (2*r*Tcn1(ny+1))+(2*(1-r)*Tcn1(ny))+(r*Tcn1(ny-1)) ;

    for j = 2:ny-2
        b(j) = (r*(Tcn1(j+2)))+(2*(1-r)*Tcn1(j+1))+(r*Tcn1(j)) ;
    end

    [X] = Tri_Diagonal(A,ny-1,b) ;

    for j = 1:ny-1
        Tcn1n(j+1) = X(j) ;
    end

    T = Tcn1n ;

end
        
           

%Boundary Condition 2 
Tft2 = T_amb*ones(ny+1,1);
Tft2(ny+1) = T_top;
Tft2n = T_amb*ones(ny+1,1);
Tft2n(ny+1) = T_top;

Tdu2 = T_amb*ones(ny+1,1);
Tdu2(ny+1) = T_top;
Tdu2n = T_amb*ones(ny+1,1);
Tdu2n(ny+1) = T_top;

Tim2 = T_amb*ones(ny+1,1);
Tim2(ny+1) = T_top;
Tim2n = T_amb*ones(ny+1,1);
Tim2n(ny+1) = T_top;

Tcn2 = T_amb*ones(ny+1,1);
Tcn2(ny+1) = T_top;
Tcn2n = T_amb*ones(ny+1,1);
Tcn2n(ny+1) = T_top;
  
 %-----FTCS Scheme-----%
        for i = 1:nt+1
            Tft2n(1) = (Tft2(2)+(300*s))/(1+s) ;
            for j = 2:ny
                Tft2n(j) = Tft2(j) + (r*(Tft2(j+1)-(2*Tft2(j))+Tft2(j-1))) ;
            end
            Tft2 = Tft2n ;
        end
        
%-----Du-Fort Frenkel Scheme-----%
        
        T_old = 300*ones(ny+1,1) ;
        T_old(ny+1) = 500 ;
                
        for i = 1:nt+1
            Tdu2n(1) = (Tdu2(2)+(300*s))/(1+s) ;
            for j = 2:ny
                Tdu2n(j) = ( T_old(j) + (2*r*(Tdu2(j+1)-T_old(j)+Tdu2(j-1))) )/(1+(2*r)) ;
            end
            T_old = Tdu2 ;
            Tdu = Tdu2n ;
        end
        
 %-----Implicit Scheme-----%
        
        A = zeros(ny-1) ;
        b = zeros(ny-1,1) ;
        
        for i = 1:nt+1
            
            Tim2n(1) = (Tim2(2)+(300*s))/(1+s) ;
            
            for j = 1:ny-1
                A(j,j) = -(1+(2*r)) ;
            end
            
            for j = 1:ny-2
                A(j,j+1) = r ;
            end
            
            for j = 2:ny-1
                A(j,j-1) = r ;
            end
            
            b(1) = -Tim2(2)-(Tim2(1)*r) ;
            b(ny-1) = -Tim2(ny)-(Tim2(ny+1)*r) ;
            
            for j = 2:ny-2
                b(j) = -Tim2(j+1) ;
            end
            
            [X] = Tri_Diagonal(A,ny-1,b) ;
            
            for j = 1:ny-1
                Tim2n(j+1) = X(j) ;
            end
            
            Tim2 = Tim2n ;
            
        end
        
 %-----Crank-Nicholson Scheme-----%
        
        A = zeros(ny-1) ;
        b = zeros(ny-1,1) ;
        
        for i = 1:nt+1
            
            Tcn2n(1) = (Tcn2(2)+(300*s))/(1+s) ;
            
            for j = 1:ny-1
                A(j,j) = 2*(1+r) ;
            end
            
            for j = 1:ny-2
                A(j,j+1) = -r ;
            end
            
            for j = 2:ny-1
                A(j,j-1) = -r ;
            end
            
            b(1) = (r*Tcn2(3))+(2*(1-r)*Tcn2(2)+(2*r*Tcn2(1))) ;
            b(ny-1) = (2*r*Tcn2(ny+1))+(2*(1-r)*Tcn2(ny))+(r*Tcn2(ny-1)) ;
            
            for j = 2:ny-2
                b(j) = (r*(Tcn2(j+2)))+(2*(1-r)*Tcn2(j+1))+(r*Tcn2(j)) ;
            end
            
            [X] = Tri_Diagonal(A,ny-1,b) ;
            
            for j = 1:ny-1
                Tcn2n(j+1) = X(j) ;
            end
            
            Tcn2 = Tcn2n ;
            
        end
            