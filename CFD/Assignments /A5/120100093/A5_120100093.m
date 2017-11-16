clear all
clc

a = 0;
b = 1;
n1 = 10;
n2 = 100;
h1 = (b-a)/n1;
h2 = (b-a)/n2;
x1 = 0:h1:1;
x2 = 0:h2:1;
y1 = zeros(1,n1+1);
y2 = zeros(1,n2+1);

for i=1:n1+1
y1(1,i) = f(x1(1,i));
end

for i=1:n2+1
y2(1,i) = f(x2(1,i));
end

% Trapezoidal Rule 

Tn1 = h1*(y1(1,1)+y1(1,n1+1))/2;
for i=2:n1
    Tn1 = Tn1 + h1*(y1(1,i));
end

Tn2 = h2*(y2(1,1)+y2(1,n2+1))/2;
for i=2:n2
    Tn2 = Tn2 + h2*(y2(1,i));
end

% Modified Trapezoidal Rule 

Tnm1 = Tn1 - ((h1*h1)/12)*(((f(x1(1,n1+1))-f(x1(1,n1)))/h1)-((f(x1(1,2))-f(x1(1,1)))/h1));
Tnm2 = Tn2 - ((h2*h2)/12)*(((f(x2(1,n2+1))-f(x2(1,n2)))/h2)-((f(x2(1,2))-f(x2(1,1)))/h2));

% Simsons Rule 

s1 = 0;
for i=2:2:n1
    s1 = s1 + y1(1,i);
end

s2 = 0;
for i=3:2:n1-1
    s2 = s2 + y1(1,i);
end

Sn1 = (h1/3)*(y1(1,1)+4*s1+2*s2+y1(1,n1+1));

s1 = 0;
for i=2:2:n2
    s1 = s1 + y2(1,i);
end

s2 = 0;
for i=3:2:n2-1
    s2 = s2 + y2(1,i);
end

Sn2 = (h2/3)*(y2(1,1)+4*s1+2*s2+y2(1,n2+1));

% Gauss Quadrature 
 
g1=(b-a)*f((b-a)/2*0+(a+b)/2); %n=1

%n=2
c=sqrt(1/3);
g2=(b-a)/2*(f((b-a)/2*c+(a+b)/2)+f((b-a)/2*(-c)+(a+b)/2)); 

%n=4
c1=sqrt((3/7)-((2/7)*(6/5)^(0.5)));
c2=sqrt((3/7)+((2/7)*(6/5)^(0.5)));
w1=(18+30^(0.5))/36;
w2=(18-30^(0.5))/36;
g3=(b-a)/2*(w1*(f((b-a)/2*c1+(a+b)/2)+f((b-a)/2*(-c1)+(a+b)/2))+w2*(f((b-a)/2*c2+(a+b)/2)+f((b-a)/2*(-c2)+(a+b)/2)));

% derivatives

c1 = 0.1;
c2 = 0.5;
yd1a = der(c1,h1);
yd1b = der(c2,h1);
yd2a = der(c1,h2);
yd2b = der(c2,h2);









    
