clc
clear all

dt = input('Enter the time-step: ');
dx = input('Enter the space-step: ');

tic

dy = dx ;
dz = dx ;

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

T1 = T ;
T2 = T ;
Tnew = T ;

rx = alpha*dt/(dx^2) ;
ry = alpha*dt/(dy^2) ;
rz = alpha*dt/(dz^2) ;

err = 15;
k = 0;

while err > 1
    
    for j = 2:ny
        for k = 2:nz
            
            A_1 = zeros(nx-1);
            b_1 = zeros(nx-1,1);
            
            for l = 2:nx-2
                A_1(l,l-1) = -rx/3 ;
                A_1(l,l) = 1+(2*rx/3) ;
                A_1(l,l+1) = -rx/3 ;
                b_1(l) = T(l+1,j,k) + (ry*(T(l+1,j+1,k)-(2*T(l+1,j,k))+T(l+1,j-1,k))/3) + (rz*(T(l+1,j,k+1)-(2*T(l+1,j,k))+T(l+1,j,k-1))/3) ;
            end
            
            A_1(1,1) = 1+(2*rx/3) ;
            A_1(1,2) = -rx/3 ;
            b_1(1) = (rx*T(1,j,k)/3) + T(2,j,k) + (ry*(T(2,j+1,k)-(2*T(2,j,k))+T(2,j-1,k))/3) + (rz*(T(2,j,k+1)-(2*T(2,j,k))+T(2,j,k-1))/3) ;
            
            A_1(nx-1,nx-1) = 1+(2*rx/3) ;
            A_1(nx-1,nx-2) = -rx/3 ;
            b_1(nx-1) = (rx*T(nx+1,j,k)/3) + T(nx,j,k) + (ry*(T(nx,j+1,k)-(2*T(nx,j,k))+T(nx,j-1,k))/3) + (rz*(T(nx,j,k+1)-(2*T(nx,j,k))+T(nx,j,k-1))/3) ;
            
            S_1 = Tri_Diagonal_Matrix_Solution(A_1,nx-1,b_1) ;
            
            for i = 2:nx
                T1(i,j,k) = S_1(i-1) ;
            end
        end
    end
    
    for i = 2:nx
        for k = 2:nz
            
            A_2 = zeros(ny-1) ;
            b_2 = zeros(ny-1,1) ;
            
            for l = 2:ny-2
                A_2(l,l-1) = -ry/3 ;
                A_2(l,l) = 1+(2*ry/3) ;
                A_2(l,l+1) = -ry/3 ;
                b_2(l) = T1(i,l+1,k) + (rx*(T1(i+1,l+1,k)-(2*T1(i,l+1,k))+T1(i-1,l+1,k))/3) + (rz*(T1(i,l+1,k+1)-(2*T1(i,l+1,k))+T1(i,l+1,k-1))/3) ;
            end
            
            A_2(1,1) = 1+(2*ry/3) ;
            A_2(1,2) = -ry/3 ;
            b_2(1) = (ry*T1(i,1,k)/3) + T1(i,2,k) + (rx*(T1(i+1,2,k)-(2*T1(i,2,k))+T1(i-1,2,k))/3) + (rz*(T1(i,2,k+1)-(2*T1(i,2,k))+T1(i,2,k-1))/3) ;
            
            A_2(ny-1,ny-1) = 1+(2*ry/3) ;
            A_2(ny-1,ny-2) = -ry/3 ;
            b_2(ny-1) = (ry*T1(i,ny+1,k)/3) + T1(i,ny,k) + (rx*(T1(i+1,ny,k)-(2*T1(i,ny,k))+T1(i-1,ny,k))/3) + (rz*(T1(i,ny,k+1)-(2*T1(i,ny,k))+T1(i,ny,k-1))/3) ;
            
            S_2 = Tri_Diagonal_Matrix_Solution(A_2,ny-1,b_2) ;
            
            for j = 2:ny
                T2(i,j,k) = S_2(j-1) ;
            end  
        end
    end
    
    for i = 2:nx
        for j = 2:ny
            
            A = zeros(nz-1) ;
            b = zeros(nz-1,1) ;
            
            for l = 2:nz-2
                A(l,l-1) = -rz/3 ;
                A(l,l) = 1+(2*rz/3) ;
                A(l,l+1) = -rz/3 ;
                b(l) = T2(i,j,l+1) + (rx*(T2(i+1,j,l+1)-(2*T2(i,j,l+1))+T2(i-1,j,l+1))/3) + (ry*(T2(i,j+1,l+1)-(2*T2(i,j,l+1))+T2(i,j-1,l+1))/3) ;
            end
            
            A(1,1) = 1+(2*rz/3) ;
            A(1,2) = -rz/3 ;
            b(1) = (rz*T2(i,j,1)/3) + T2(i,j,2) + (rx*(T2(i+1,j,2)-(2*T2(i,j,2))+T2(i-1,j,2))/3) + (ry*(T2(i,j+1,2)-(2*T2(i,j,2))+T2(i,j-1,2))/3) ;
            
            A(nz-1,nz-1) = 1+(2*rz/3) ;
            A(nz-1,nz-2) = -rz/3 ;
            b(nz-1) = (rz*T2(i,j,nz+1)/3) + T2(i,j,nz) + (rx*(T2(i+1,j,nz)-(2*T2(i,j,nz))+T2(i-1,j,nz))/3) + (ry*(T2(i,j+1,nz)-(2*T2(i,j,nz))+T2(i,j-1,nz))/3) ;
            
            S = Tri_Diagonal_Matrix_Solution(A,nz-1,b) ;
            
            for k = 2:nz
                Tnew(i,j,k) = S(k-1) ;
            end
        end
    end
    
    k = k+1;
    T_measured(k,1) = Tnew((0.5/dx)+1,fix(0.25/dy)+1,fix(0.25/dz)+1) ;
    T_measured(k,2) = Tnew((1/dx)+1,(0.5/dy)+1,(0.5/dz)+1) ;
    T_measured(k,3) = Tnew((1.5/dx)+1,fix(0.75/dy)+1,fix(0.75/dz)+1) ;
    T_measured(k,4) = Tnew((0.5/dx)+1,(0.5/dy)+1,(0.5/dz)+1) ;
    T_measured(k,5) = Tnew((1/dx)+1,fix(0.75/dy)+1,fix(0.75/dz)+1) ;
    
    Time(k,1) = (k-1)*dt ;
    err = max(max(max(abs(Tnew-T))))
    T = Tnew ;
end

timeTaken = toc;

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

figure
plot(Time,T_measured)
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