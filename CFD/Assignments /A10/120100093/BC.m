function[T]=BC(Ti)
[N,M]=size(Ti);
T=Ti;
for i=1:N
    T(i,1)=100;
    T(i,M)=1000;
end

for j=1:M
    T(1,j)=0;
    T(N,j)=500;
end
end