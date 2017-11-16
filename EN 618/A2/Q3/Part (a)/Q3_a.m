clc
clear all
format long

delta_t = input('Enter the value of time-step: ') ;
delta_x = input('Enter the value of space-step: ') ;

L = 5 ;
Time = 8;
alpha = 10^(-4); %input('Enter the value of alpha: ') ;
x0 = 0.5 ;
sigma = 0.01 ;
u = 0.2 ;

X = 0:delta_x:L ;

n_t = Time/delta_t ;
n_x = L/delta_x ;

r = alpha*delta_t/(delta_x^2) ;
q = u*delta_t/delta_x ;

T = zeros(n_x+1,1) ;
T_new = zeros(n_x+1,1) ;

for i = 2:n_x
    T(i) = exp(-((X(i)-x0)^2)/sigma) ;
end

%-----Exact Function-----%

T_exact = zeros(n_x+1,1) ;

for i=1:n_x+1
    T_exact(i) = T_exact_function(X(i),Time,alpha,sigma,x0,u);
end

%-----Explicit Euler-----%

T_e = T ;
for i = 1:n_t
    for j = 2:n_x
        T_new(j) = T_e(j) + (r*(T_e(j+1)-(2*T_e(j))+T_e(j-1))) + (q*(T_e(j+1)-T_e(j-1))/2) ;
    end
    T_e = T_new ;
end

%-----Implicit Euler-----%

T_i = T ;

A_i = zeros(n_x-1) ;
b_i = zeros(n_x-1,1) ;

for k_i=1:n_t
    
    for i = 2:n_x-2
        A_i(i,i-1) = -r-(q/2) ;
        A_i(i,i) = (2*r)+1 ;
        A_i(i,i+1) = -r+(q/2) ;
        b_i(i) = T_i(i+1) ;
    end
    
    A_i(1,1) = (2*r)+1 ;
    A_i(1,2) = -r+(q/2) ;
    b_i(1) = T_i(2)+((r+(q/2))*T(1)) ;

    A_i(n_x-1,n_x-1) = (2*r)+1 ;
    A_i(n_x-1,n_x-2) = -r-(q/2) ;
    b_i(n_x-1) = T_i(n_x)+((-r+(q/2))*T(n_x+1)) ;

    [A_i,S_i,b_i]=Tri_Diagonal_Matrix_Solution(A_i,n_x-1,b_i);
    
    for i=1:n_x-1
        T_i(i+1) = S_i(i) ;
    end
    k_i
end

%-----Trapezoidal-----%

T_t = T ;

A_t = zeros(n_x-1) ;
b_t = zeros(n_x-1,1) ;

for k_t=1:n_t
    
    for i = 2:n_x-2
        A_t(i,i-1) = -(r/2)-(q/4) ;
        A_t(i,i) = r+1 ;
        A_t(i,i+1) = -(r/2)+(q/4) ;
    end
    
    A_t(1,1) = r+1 ;
    A_t(1,2) = -(r/2)+(q/4) ;
    A_t(n_x-1,n_x-1) = r+1 ;
    A_t(n_x-1,n_x-2) = -(r/2)-(q/4) ;

    for i = 2:n_x-2
        b_t(i) = (r*(T_t(i+2)-(2*T_t(i+1))+T_t(i))/2) + T_t(i+1) - (q*(T_t(i+1)-T_t(i-1))/4) ;
    end
    
    b_t(1) = (r*(T_t(3)-(2*T_t(2))+T_t(1))/2) + T_t(2) + (((r/2)+(q/4))*T_t(1)) ;
    b_t(n_x-1) = (r*(T_t(n_x+1)-(2*T_t(n_x))+T_t(n_x-1))/2) + T_t(n_x) + (((r/2)-(q/4))*T_t(n_x+1)) ;
    
    [A_t,S_t,b_t]=Tri_Diagonal_Matrix_Solution(A_t,n_x-1,b_t);
    
    for i=1:n_x-1
        T_t(i+1) = S_t(i) ;
    end
    k_t
end

%-----Runge Kutta 4th Order Scheme-----%

T_r = T ;
k1 = zeros (n_x-1,1) ;
k2 = zeros (n_x-1,1) ;
k3 = zeros (n_x-1,1) ;
k4 = zeros (n_x-1,1) ;

for k_r = 1:n_t
    
    for i=1:n_x-1
        k1(i) = ( alpha * ( T_r(i+2) - (2*T_r(i+1)) + T_r(i) ) / (delta_x^2) ) - ( u * ( T_r(i+2) - T_r(i) ) / (2*delta_x )) ;
    end
    
    for i=1:n_x-1
        if i==1
            k2(i) = alpha* ( (T_r(i+2)+(delta_t*k1(i+1)/2)) - (2*(T_r(i+1)+(delta_t*k1(i)/2))) + T_r(i) ) / (delta_x^2) - ( u * ( (T_r(i+2)+(delta_t*k1(i+1)/2)) - T_r(i) ) / (2*delta_x) );
        elseif i==n_x-1
            k2(i) = alpha* ( T_r(i+2) - (2*(T_r(i+1)+(delta_t*k1(i)/2))) + (T_r(i)+(delta_t*k1(i-1)/2)) ) / (delta_x^2) - ( u * ( T_r(i+2) - (T_r(i)+(delta_t*k1(i-1)/2)) ) / (2*delta_x) );
        else
            k2(i) = alpha* ( (T_r(i+2)+(delta_t*k1(i+1)/2)) - (2*(T_r(i+1)+(delta_t*k1(i)/2))) + (T_r(i)+(delta_t*k1(i-1)/2)) ) / (delta_x^2) - ( u * ( (T_r(i+2)+(delta_t*k1(i+1)/2)) - (T_r(i)+(delta_t*k1(i-1)/2)) ) / (2*delta_x) ) ;
        end
    end
    
    for i=1:n_x-1
        if i==1
            k3(i) = alpha* ( (T_r(i+2)+(delta_t*k2(i+1)/2)) - (2*(T_r(i+1)+(delta_t*k2(i)/2))) + T_r(i) ) / (delta_x^2) - ( u * ( (T_r(i+2)+(delta_t*k2(i+1)/2)) - T_r(i) ) / (2*delta_x) );
        elseif i==n_x-1
            k3(i) = alpha* ( T_r(i+2) - (2*(T_r(i+1)+(delta_t*k2(i)/2))) + (T_r(i)+(delta_t*k2(i-1)/2)) ) / (delta_x^2) - ( u * ( T_r(i+2) - (T_r(i)+(delta_t*k2(i-1)/2)) ) / (2*delta_x) );
        else
            k3(i) = alpha* ( (T_r(i+2)+(delta_t*k2(i+1)/2)) - (2*(T_r(i+1)+(delta_t*k2(i)/2))) + (T_r(i)+(delta_t*k2(i-1)/2)) ) / (delta_x^2) - ( u * ( (T_r(i+2)+(delta_t*k2(i+1)/2)) - (T_r(i)+(delta_t*k2(i-1)/2)) ) / (2*delta_x) );
        end
    end
    
    for i=1:n_x-1
        if i==1
            k4(i) = alpha* ( (T_r(i+2)+(delta_t*k3(i+1))) - (2*(T_r(i+1)+(delta_t*k3(i)))) + T_r(i) ) / (delta_x^2) - ( u * ( (T_r(i+2)+(delta_t*k3(i+1))) - T_r(i) ) / (2*delta_x) );
        elseif i==n_x-1
            k4(i) = alpha* ( T_r(i+2) - (2*(T_r(i+1)+(delta_t*k3(i)))) + (T_r(i)+(delta_t*k3(i-1))) ) / (delta_x^2) - ( u * ( T_r(i+2) - (T_r(i)+(delta_t*k3(i-1))) ) / (2*delta_x) );
        else
            k4(i) = alpha* ( (T_r(i+2)+(delta_t*k3(i+1))) - (2*(T_r(i+1)+(delta_t*k3(i)))) + (T_r(i)+(delta_t*k3(i-1))) ) / (delta_x^2) - ( u * ( (T_r(i+2)+(delta_t*k3(i+1))) - (T_r(i)+(delta_t*k3(i-1))) ) / (2*delta_x) );
        end
    end
    
    for i=1:n_x-1
        T_new(i+1) = T_r(i+1) + (delta_t*(k1(i)+(2*(k2(i)+k3(i)))+k4(i))/6) ;
    end
    
    T_r = T_new ;
    
    k_r
    
end