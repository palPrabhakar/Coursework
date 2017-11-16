clear all
clc

h=[0.1 0.01 0.001];
N=[0 1 10];
x_initial=0;
y_initial=0;
x_max=5;

yf1_1=f1(x_initial,y_initial,x_max,h(1)); %Exact function value
yf1_2=f1(x_initial,y_initial,x_max,h(2)); %Exact function value
yf1_3=f1(x_initial,y_initial,x_max,h(3)); %Exact function value
yf21_1=f2(x_initial,y_initial,x_max,h(1),-1); %Exact function value
yf21_2=f2(x_initial,y_initial,x_max,h(2),-1); %Exact function value
yf21_3=f2(x_initial,y_initial,x_max,h(3),-1); %Exact function value
yf22_1=f2(x_initial,y_initial,x_max,h(1),-50); %Exact function value
yf22_2=f2(x_initial,y_initial,x_max,h(2),-50); %Exact function value
yf22_3=f2(x_initial,y_initial,x_max,h(3),-50); %Exact function value

for j=1:3
    n(j)=(x_max-x_initial)/h(j); 
end

x_1(1)=x_initial;
x_2(1)=x_initial;
x_3(1)=x_initial;

for i=1:n(1)
    x_1(i+1)=x_1(i)+h(1);
end

for i=1:n(2)
    x_2(i+1)=x_2(i)+h(2);
end

for i=1:n(3)
    x_3(i+1)=x_3(i)+h(3);
end

[yfe1_1,yfe21_1,yfe22_1]=forward_euler(y_initial,x_1,h(1),n(1),-1,-50);
[yfe1_2,yfe21_2,yfe22_2]=forward_euler(y_initial,x_2,h(2),n(2),-1,-50);
[yfe1_3,yfe21_3,yfe22_3]=forward_euler(y_initial,x_3,h(3),n(3),-1,-50); % function  1 and function 2 with lambda=-1 and lambda=-50 for h=0.001

[ybe1_1_n1,ybe21_1_n1,ybe22_1_n1]=backward_euler(y_initial,yfe1_1,yfe21_1,yfe22_1,x_1,h(1),n(1),0,-1,-50);
[ybe1_2_n1,ybe21_2_n1,ybe22_2_n1]=backward_euler(y_initial,yfe1_2,yfe21_2,yfe22_2,x_2,h(2),n(2),0,-1,-50);
[ybe1_3_n1,ybe21_3_n1,ybe22_3_n1]=backward_euler(y_initial,yfe1_3,yfe21_3,yfe22_3,x_3,h(3),n(3),0,-1,-50);

[ybe1_1_n2,ybe21_1_n2,ybe22_1_n2]=backward_euler(y_initial,yfe1_1,yfe21_1,yfe22_1,x_1,h(1),n(1),1,-1,-50);
[ybe1_2_n2,ybe21_2_n2,ybe22_2_n2]=backward_euler(y_initial,yfe1_2,yfe21_2,yfe22_2,x_2,h(2),n(2),1,-1,-50);
[ybe1_3_n2,ybe21_3_n2,ybe22_3_n2]=backward_euler(y_initial,yfe1_3,yfe21_3,yfe22_3,x_3,h(3),n(3),1,-1,-50);

[ybe1_1_n3,ybe21_1_n3,ybe22_1_n3]=backward_euler(y_initial,yfe1_1,yfe21_1,yfe22_1,x_1,h(1),n(1),10,-1,-50);
[ybe1_2_n3,ybe21_2_n3,ybe22_2_n3]=backward_euler(y_initial,yfe1_2,yfe21_2,yfe22_2,x_2,h(2),n(2),10,-1,-50);
[ybe1_3_n3,ybe21_3_n3,ybe22_3_n3]=backward_euler(y_initial,yfe1_3,yfe21_3,yfe22_3,x_3,h(3),n(3),10,-1,-50);


[yt1_1_n1,yt21_1_n1,yt22_1_n1]=trapezoidal(y_initial,yfe1_1,yfe21_1,yfe22_1,x_1,h(1),n(1),0,-1,-50);
[yt1_2_n1,yt21_2_n1,yt22_2_n1]=trapezoidal(y_initial,yfe1_2,yfe21_2,yfe22_2,x_2,h(2),n(2),0,-1,-50);
[yt1_3_n1,yt21_3_n1,yt22_3_n1]=trapezoidal(y_initial,yfe1_3,yfe21_3,yfe22_3,x_3,h(3),n(3),0,-1,-50);

[yt1_1_n2,yt21_1_n2,yt22_1_n2]=trapezoidal(y_initial,yfe1_1,yfe21_1,yfe22_1,x_1,h(1),n(1),1,-1,-50);
[yt1_2_n2,yt21_2_n2,yt22_2_n2]=trapezoidal(y_initial,yfe1_2,yfe21_2,yfe22_2,x_2,h(2),n(2),1,-1,-50);
[yt1_3_n2,yt21_3_n2,yt22_3_n2]=trapezoidal(y_initial,yfe1_3,yfe21_3,yfe22_3,x_3,h(3),n(3),1,-1,-50);

[yt1_1_n3,yt21_1_n3,yt22_1_n3]=trapezoidal(y_initial,yfe1_1,yfe21_1,yfe22_1,x_1,h(1),n(1),10,-1,-50);
[yt1_2_n3,yt21_2_n3,yt22_2_n3]=trapezoidal(y_initial,yfe1_2,yfe21_2,yfe22_2,x_2,h(2),n(2),10,-1,-50);
[yt1_3_n3,yt21_3_n3,yt22_3_n3]=trapezoidal(y_initial,yfe1_3,yfe21_3,yfe22_3,x_3,h(3),n(3),10,-1,-50);


%{
plot(x_3,yf1_3,'b');
hold on 
plot(x_3,yfe1_3,'r');
hold on 
plot(x_3,ybe1_3_n1,'g');
plot(x_3,yt1_3_n1,'y');
title('Plot of first function for h=0.0001 and N=0');
legend('Exact Solution','Forward Euler','Backward Euler','Trapezoidal');
%}

%{
plot(x_3,yf1_3,'b');
hold on 
plot(x_3,yfe1_3,'r');
hold on 
plot(x_3,ybe1_3_n2,'g');
plot(x_3,yt1_3_n2,'y');
title('Plot of first function for h=0.0001 and N=1');
legend('Exact Solution','Forward Euler','Backward Euler','Trapezoidal');
%}

%{
plot(x_3,yf1_3,'b');
hold on 
plot(x_3,yfe1_3,'r');
hold on 
plot(x_3,ybe1_3_n3,'g');
plot(x_3,yt1_3_n3,'y');
title('Plot of first function for h=0.0001 and N=10');
legend('Exact Solution','Forward Euler','Backward Euler','Trapezoidal');
%}

%{
plot(x_3,yf1_3,'b')
hold on 
plot(x_1,yfe1_1,'r') 
plot(x_2,yfe1_2,'g')
plot(x_3,yfe1_3,'y')
title('Plot of first function using forward euler method')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf1_3,'b')
hold on 
plot(x_1,ybe1_1_n1,'r') 
plot(x_2,ybe1_2_n1,'g')
plot(x_3,ybe1_3_n1,'y')
title('Plot of first function using backward euler method and N=0')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}


%{
plot(x_3,yf1_3,'b')
hold on 
plot(x_1,ybe1_1_n2,'r') 
plot(x_2,ybe1_2_n2,'g')
plot(x_3,ybe1_3_n2,'y')
title('Plot of first function using backward euler method and N=1')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf1_3,'b')
hold on 
plot(x_1,ybe1_1_n3,'r') 
plot(x_2,ybe1_2_n3,'g')
plot(x_3,ybe1_3_n3,'y')
title('Plot of first function using backward euler method and N=10')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf1_3,'b')
hold on 
plot(x_1,yt1_1_n1,'r') 
plot(x_2,yt1_2_n1,'g')
plot(x_3,yt1_3_n1,'y')
title('Plot of first function using trapezoidal method and N=0')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf1_3,'b')
hold on 
plot(x_1,yt1_1_n2,'r') 
plot(x_2,yt1_2_n2,'g')
plot(x_3,yt1_3_n2,'y')
title('Plot of first function using trapezoidal method and N=1')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf1_3,'b')
hold on 
plot(x_1,yt1_1_n3,'r') 
plot(x_2,yt1_2_n3,'g')
plot(x_3,yt1_3_n3,'y')
title('Plot of first function using trapezoidal method and N=10')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf21_3,'b');
hold on 
plot(x_3,yfe21_3,'r');
hold on 
plot(x_3,ybe21_3_n1,'g');
plot(x_3,yt21_3_n1,'y');
title('Plot of second function for h=0.0001 and N=0 and lambda=-1');
legend('Exact Solution','Forward Euler','Backward Euler','Trapezoidal');
%}


%{
plot(x_3,yf22_3,'b');
hold on 
plot(x_3,yfe22_3,'r');
hold on 
plot(x_3,ybe22_3_n1,'g');
plot(x_3,yt22_3_n1,'y');
title('Plot of second function for h=0.0001 and N=0 and lambda=-50');
legend('Exact Solution','Forward Euler','Backward Euler','Trapezoidal');
%}


%{
plot(x_3,yf21_3,'b');
hold on 
plot(x_3,yfe21_3,'r');
hold on 
plot(x_3,ybe21_3_n1,'g');
plot(x_3,yt21_3_n1,'y');
title('Plot of second function for h=0.0001 and N=0 and lambda=-1');
legend('Exact Solution','Forward Euler','Backward Euler','Trapezoidal');
%}

%{
plot(x_3,yf22_3,'b');
hold on 
plot(x_3,yfe22_3,'r');
hold on 
plot(x_3,ybe22_3_n1,'g');
plot(x_3,yt22_3_n1,'y');
title('Plot of second function for h=0.0001 and N=0 and lambda=-50');
legend('Exact Solution','Forward Euler','Backward Euler','Trapezoidal');
%}

%{
plot(x_3,yf21_3,'b');
hold on 
plot(x_3,yfe21_3,'r');
hold on 
plot(x_3,ybe21_3_n2,'g');
plot(x_3,yt21_3_n2,'y');
title('Plot of second function for h=0.0001 and N=1 and lambda=-1');
legend('Exact Solution','Forward Euler','Backward Euler','Trapezoidal');
%}


%{
plot(x_3,yf21_3,'b');
hold on 
plot(x_3,yfe21_3,'r');
hold on 
plot(x_3,ybe21_3_n3,'g');
plot(x_3,yt21_3_n3,'y');
title('Plot of second function for h=0.0001 and N=10 and lambda=-1');
legend('Exact Solution','Forward Euler','Backward Euler','Trapezoidal');
%}

%{
plot(x_3,yf22_3,'b');
hold on 
plot(x_3,yfe22_3,'r');
hold on 
plot(x_3,ybe22_3_n2,'g');
plot(x_3,yt22_3_n2,'y');
title('Plot of second function for h=0.0001 and N=1 and lambda=-50');
legend('Exact Solution','Forward Euler','Backward Euler','Trapezoidal');
%}


%{
plot(x_3,yf22_3,'b');
hold on 
plot(x_3,yfe22_3,'r');
hold on 
plot(x_3,ybe22_3_n3,'g');
plot(x_3,yt22_3_n3,'y');
title('Plot of second function for h=0.0001 and N=10 and lambda=-50');
legend('Exact Solution','Forward Euler','Backward Euler','Trapezoidal');
%}

%{
plot(x_3,yf21_3,'b')
hold on 
plot(x_1,yfe21_1,'r') 
plot(x_2,yfe21_2,'g')
plot(x_3,yfe21_3,'y')
title('Plot of second function using forward euler method and lambda=-1')
legend('Exact Solution','h=0.1','h=0.01','h=0.001');
%}

%{
plot(x_3,yf22_3,'b')
hold on 
plot(x_1,yfe22_1,'r') 
plot(x_2,yfe22_2,'g')
plot(x_3,yfe22_3,'y')
title('Plot of second function using forward euler method and lambda=-50')
legend('Exact Solution','h=0.1','h=0.01','h=0.001');
%}

%{
plot(x_3,yf21_3,'b')
hold on 
plot(x_1,ybe21_1_n1,'r') 
plot(x_2,ybe21_2_n1,'g')
plot(x_3,ybe21_3_n1,'y')
title('Plot of second function using backward euler method and N=0 and lambda=-1')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}


%{
plot(x_3,yf21_3,'b')
hold on 
plot(x_1,ybe21_1_n2,'r') 
plot(x_2,ybe21_2_n2,'g')
plot(x_3,ybe21_3_n2,'y')
title('Plot of second function using backward euler method and N=1 and lambda=-1')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf21_3,'b')
hold on 
plot(x_1,ybe21_1_n3,'r') 
plot(x_2,ybe21_2_n3,'g')
plot(x_3,ybe21_3_n3,'y')
title('Plot of second function using backward euler method and N=10 and lambda=-1')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf22_3,'b')
hold on 
plot(x_1,ybe22_1_n1,'r') 
plot(x_2,ybe22_2_n1,'g')
plot(x_3,ybe22_3_n1,'y')
title('Plot of second function using backward euler method and N=0 and lambda=-50')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}


%{
plot(x_3,yf22_3,'b')
hold on 
plot(x_1,ybe22_1_n2,'r') 
plot(x_2,ybe22_2_n2,'g')
plot(x_3,ybe22_3_n2,'y')
title('Plot of second function using backward euler method and N=1 and lambda=-50')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf22_3,'b')
hold on 
plot(x_1,ybe22_1_n3,'r') 
plot(x_2,ybe22_2_n3,'g')
plot(x_3,ybe22_3_n3,'y')
title('Plot of second function using backward euler method and N=10 and lambda=-50')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf21_3,'b')
hold on 
plot(x_1,yt21_1_n1,'r') 
plot(x_2,yt21_2_n1,'g')
plot(x_3,yt21_3_n1,'y')
title('Plot of second function using trapezoidal method and N=0 and lambda=-1')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf21_3,'b')
hold on 
plot(x_1,yt21_1_n2,'r') 
plot(x_2,yt21_2_n2,'g')
plot(x_3,yt21_3_n2,'y')
title('Plot of seocnd function using trapezoidal method and N=1 and l=-1')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf21_3,'b')
hold on 
plot(x_1,yt21_1_n3,'r') 
plot(x_2,yt21_2_n3,'g')
plot(x_3,yt21_3_n3,'y')
title('Plot of 2 function using trapezoidal method and N=10 and l=-1')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf22_3,'b')
hold on 
plot(x_1,yt22_1_n1,'r') 
plot(x_2,yt22_2_n1,'g')
plot(x_3,yt22_3_n1,'y')
title('Plot of second function using trapezoidal method and N=0 and lambda=-50')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf22_3,'b')
hold on 
plot(x_1,yt22_1_n2,'r') 
plot(x_2,yt22_2_n2,'g')
plot(x_3,yt22_3_n2,'y')
title('Plot of seocnd function using trapezoidal method and N=1 and l=-50')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}

%{
plot(x_3,yf22_3,'b')
hold on 
plot(x_1,yt22_1_n3,'r') 
plot(x_2,yt22_2_n3,'g')
plot(x_3,yt22_3_n3,'y')
title('Plot of 2 function using trapezoidal method and N=10 and l=-50')
legend('Exact Solution','h=0.1','h=0.01','h=0.001')
%}
