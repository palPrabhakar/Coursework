function [err] = ric_err(y,yrk2h,yrk2h2,yrk3h,yrk3h2,yrk4h,yrk4h2,yab4h,yab4h2,h,n1,n2)

k1=2;
k2=3;
k3=4;
t1=n1/5;
t2=n2/5;

for i=1:5
    for j=1:4
        err(i,j)=y(1+t1*i)-((((2^k1)*yrk2h2(1+t2*i))-yrk2h(1+t1*i))/((2^k1)-1));
        err(i,j)=y(1+t1*i)-((((2^k2)*yrk3h2(1+t2*i))-yrk3h(1+t1*i))/((2^k2)-1));
        err(i,j)=y(1+t1*i)-((((2^k3)*yrk4h2(1+t2*i))-yrk4h(1+t1*i))/((2^k3)-1));
        err(i,j)=y(1+t1*i)-((((2^k3)*yab4h2(1+t2*i))-yab4h(1+t1*i))/((2^k3)-1));
    end
end
