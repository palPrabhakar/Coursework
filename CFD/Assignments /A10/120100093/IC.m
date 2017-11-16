function[T]=IC(Ti,I,dx,dy)
[N,M]=size(Ti);
T=Ti;
if I=='A'
    for i=2:N-1
        for j=2:M-1
            T(i,j)=0;
        end
    end
end

if I=='B'
    for i=2:N-1
        for j=2:M-1
            T(i,j)=200;
        end
    end
end

end