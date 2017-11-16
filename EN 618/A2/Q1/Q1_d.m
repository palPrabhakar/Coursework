clc
clear all

n = input('Enter the number of grid points: ') ;
method = input('Enter \n1 for Gauss-Siedel Method \n2 for Gaussian Elimination\n ');

T = zeros (n,1) ;

delta_x = 1/(n-1) ;
x = zeros (n,1) ;
for i=1:n
    x(i) = (i-1)*delta_x ;
end

T_exact = zeros(n,1) ;
for i=1:n
    T_exact(i) = (1-(x(i)^2))/2 ;
end

switch method
    
    case 1
        %-----Gauss-Seidel Method-----%
        
        T_new = zeros (n,1) ;
        k=0;
        err=1;
        
        while err > 10^(-6)
            
            k=k+1;
            
            for i=2:n-1
                T_new(i) = ((delta_x^2) + T_new(i-1) + T(i+1)) / 2 ;
            end
            T_new(1) = T_new(2);
            err = abs(max( T_new - T ))
            T = T_new ;
            
            if k>2000
                break
            end
        end
    
    case 2
        %-----Gaussian Elimination-----%
        
        A = zeros(n-2) ;
        b = zeros(n-2,1) ;
        
        for i=1:n-2
            if i==1
                A(i,i) = -1;
                A(i,i+1) = 1;
                b(i) = -(delta_x^2) ;
            elseif i==n-2
                A(i,i-1) = 1;
                A(i,i) = -2;
                b(i) = -(delta_x^2)-T(i) ;
            else
                A(i,i-1) = 1;
                A(i,i) = -2;
                A(i,i+1) = 1;
                b(i) = -(delta_x^2) ;
            end
        end
        
        X =zeros(n-2,1);
        for k=1:n-2
            for i=k+1:n-2
                m(i,k)=A(i,k)/A(k,k);
                for j=k:n-2
                    A(i,j) = A(i,j) - (m(i,k)*A(k,j));
                end
                b(i) = b(i)-(m(i,k)*b(k));
            end
        end
        
        X(n-2) = b(n-2)/A(n-2,n-2);
        for k = n-3:-1:1
            sum = 0;
            for j = k+1:n-2
                sum = sum + (X(j)*A(k,j));
            end
            X(k) = (b(k)-sum)/A(k,k);
        end
        
        for i=1:n-2
            T(i+1) = X(i);
        end
        T(1)=T(2);
    otherwise fprintf('Not a valid option')
end

plot(x,T,'r--o');
 hold on 
 plot(x,T_exact,'b--*');
 legend('Numerical','Exact');
 title('Q1-d for N = 33');