clc
clear all

dt = input('Enter the time-step: ');
dx = input('Enter the space-step: ');

tic

dy = dx;
dz =dx;

nx = 2/dx ;
ny = 1/dy ;
nz = 1/dz ;

X = 0:dx:2 ;
Y = 0:dy:1 ;
Z = 0:dz:1 ;

rho = 4000 ;
Cp = 200 ;
k = 400 ;

alpha = k/(rho*Cp) ;

T = 300*ones(nx+1,ny+1,nz+1) ;

T(1,:,:) = 100 ;
T(nx+1,:,:) = 500 ;
T(:,1,:) = 300 ;
T(:,ny+1,:) = 600 ;
T(:,:,1) = 500 ;
T(:,:,nz+1) = 2000 ;

T_new = T ;

rx = alpha*dt/(dx^2) ;
ry = alpha*dt/(dy^2) ;
rz = alpha*dt/(dz^2) ;

err = 25;
k = 0;

while err>.001
    
    for i=2:nx
        for j=2:ny
            for k=2:nz
                T_new(i,j,k) = T(i,j,k) + (rx*(T(i+1,j,k)-(2*T(i,j,k))+T(i-1,j,k))) + (ry*(T(i,j+1,k)-(2*T(i,j,k))+T(i,j-1,k))) + (rz*(T(i,j,k+1)-(2*T(i,j,k))+T(i,j,k-1))) ;
            end
        end
    end
    
    k=k+1;
    T_trans(k,1) = T_new((0.5/dx)+1,fix(0.25/dy)+1,fix(0.25/dz)+1) ;
    T_trans(k,2) = T_new((1/dx)+1,(0.5/dy)+1,(0.5/dz)+1) ;
    T_trans(k,3) = T_new((1.5/dx)+1,fix(0.75/dy)+1,fix(0.75/dz)+1) ;
    T_trans(k,4) = T_new((0.5/dx)+1,(0.5/dy)+1,(0.5/dz)+1) ;
    T_trans(k,5) = T_new((1/dx)+1,fix(0.75/dy)+1,fix(0.75/dz)+1) ;
    Time(k,1) = (k-1)*dt ;
    err = max(max(max(abs(T_new-T))))
    T = T_new ;   
end

for i = (1/dx)+1
    for j = 1:ny+1
        for k = 1:nz+1
            T_X(j,k) = T(i,j,k) ;
        end
    end
end

for j = (0.5/dy)+1
    for i = 1:nx+1
        for k = 1:nz+1
            T_Y(i,k) = T(i,j,k) ;
        end
    end
end

for k = (0.5/dz)+1
    for j = 1:ny+1
        for i = 1:nx+1
            T_Z(i,j) = T(i,j,k) ;
        end
    end
end

total_time = toc;

figure
plot(Time,T_trans)
title('Temperature History')
legend('Point 1','Point 2','Point 3','Point 4','Point 5')
xlabel('Time')
ylabel('Temperature')

figure
contourf(Z,Y,T_X)
title('Temperature Profile at X=1')
xlabel('Z')
ylabel('Y')
colorbar

figure
contourf(Z,X,T_Y)
title('Temperature Profile at Y=0.5')
xlabel('Z')
ylabel('X')
colorbar

figure
contourf(Y,X,T_Z)
title('Temperature Profile at Z=0.5')
xlabel('Y')
ylabel('X')
colorbar