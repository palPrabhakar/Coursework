function [y,x] = exact(x_initial,n,h);

x(1)=x_initial;

for j=1:n
    x(j+1)=x(j)+h;
end

for i=1:n+1
    y(i)=x(i)^1.1;
end

