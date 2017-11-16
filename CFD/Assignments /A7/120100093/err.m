function [e] = err(y,y1,y2,y3,y4,h,n);

t=n/5;
j=1;
a=t+1;
b=t;
for i=a:b:n+1
        e(j,1)=y(i)-y1(i);
        e(j,2)=y(i)-y2(i);
        e(j,3)=y(i)-y3(i);
        e(j,4)=y(i)-y4(i);
        j=j+1;
end

