function [y1,y2,y3] = trapezoidal(a,y_1,y_2,y_3,x,h,n,N,l1,l2)

y1(1)=a;
y2(1)=a;
y3(1)=a;
err=1;
epsilon=1*10^-6;

for i=1:n

    if N==0
    
        if i==1
            y1(1,i+1) = y1(i) + 0.5*h*((x(i+1)*cos(x(i+1)))+(y_1(i+1)/x(i+1)));
            err=y1(i+1)-y_1(i+1);
            count=0;
            while(abs(err)>epsilon && count<100)
                temp=y1(i);
                y1(i)=y1(i)+0.5*h*((x(i+1)*cos(x(i+1)))+(y1(i+1)/x(i+1)));
                err=y1(i)-temp;
                count=count+1;
            end
        else
            y1(1,i+1) = y1(i) + 0.5*h*((x(i)*cos(x(i)))+(y1(i)/x(i))+((x(i+1)*cos(x(i+1)))+(y_1(i+1)/x(i+1))));
            err=y1(i+1)-y_1(i+1);
            count=0;
            while(abs(err)>epsilon && count<100)
                temp=y1(i);
                y1(i)=y1(i)+0.5*h*((x(i)*cos(x(i)))+(y1(i)/x(i))+((x(i+1)*cos(x(i+1)))+(y1(i+1)/x(i+1))));
                err=y1(i)-temp;
                count=count+1;
            end
        end
        
        y2(i+1)=y2(i)+0.5*h*(((l1*y2(i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i)))))+((l1*y_2(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l1*(atan(x(i+1))))));
        err=y2(i)-y_2(i);
        count=0;
        while(abs(err)>epsilon && count<100)
            temp=y2(i+1);
            y2(i+1)=y2(i)+0.5*h*(((l1*y2(i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i)))))+((l1*y2(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l1*(atan(x(i+1))))));
            err=y2(i+1)-temp;
            count=count+1;
        end
        
        y3(i+1)=y3(i)+0.5*h*(((l2*y3(i))+(1/(1+(x(i)*x(i))))-(l2*(atan(x(i)))))+((l2*y_3(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l2*(atan(x(i+1))))));
        err=y3(i)-y_3(i);
        count=0;
        while(abs(err)>epsilon && count<100)
            temp=y2(i);
            y3(i)=y3(i)+0.5*h*(((l2*y3(i))+(1/(1+(x(i)*x(i))))-(l2*(atan(x(i)))))+((l2*y3(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l2*(atan(x(i+1))))));
            err=y3(i+1)-temp;
            count=count+1;
        end
        
    end
    
    if N==1
        
        if i==1
            y1(i+1)=y1(i)+0.5*h*((x(i+1)*cos(x(i+1)))+(y_1(i+1)/x(i+1)));
            y1(i+1)=y1(i)+0.5*h*((x(i+1)*cos(x(i+1)))+(y1(i+1)/x(i+1)));
        else
            y1(i+1)=y1(i)+0.5*h*(((x(i)*cos(x(i)))+(y_1(i)/x(i)))+((x(i+1)*cos(x(i+1)))+(y_1(i+1)/x(i+1))));
            y1(i+1)=y1(i)+0.5*h*(((x(i)*cos(x(i)))+(y1(i)/x(i)))+((x(i+1)*cos(x(i+1)))+(y1(i+1)/x(i+1)))); 
        end
        
         y2(i+1)=y2(i)+0.5*h*(((l1*y2(i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i))))+((l1*y_2(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l1*(atan(x(i+1)))))));
         y2(i+1)=y2(i)+0.5*h*(((l1*y2(i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i))))+((l1*y2(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l1*(atan(x(i+1)))))));
         
         y3(i+1)=y3(i)+0.5*h*(((l2*y3(i))+(1/(1+(x(i)*x(i))))-(l2*(atan(x(i))))+((l2*y_3(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l2*(atan(x(i+1)))))));
         y3(i+1)=y3(i)+0.5*h*(((l2*y3(i))+(1/(1+(x(i)*x(i))))-(l2*(atan(x(i))))+((l2*y3(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l2*(atan(x(i+1)))))));
    
    end
    
    if N==10
    
        if i==1
            y1(i+1)=y1(i)+0.5*h*((x(i+1)*cos(x(i+1)))+(y_1(i+1)/x(i+1)));
            count=0;
            while(count<10)
                y1(i+1)=y1(i)+0.5*h*((x(i+1)*cos(x(i+1)))+(y1(i+1)/x(i+1)));
                count=count+1;
            end
        else
            y1(i+1)=y1(i)+0.5*h*(((x(i)*cos(x(i)))+(y_1(i)/x(i)))+((x(i+1)*cos(x(i+1)))+(y_1(i+1)/x(i+1))));
            count=0;
            while(count<10)
                y1(i+1)=y1(i)+0.5*h*(((x(i)*cos(x(i)))+(y1(i)/x(i)))+((x(i+1)*cos(x(i+1)))+(y1(i+1)/x(i+1)))); 
                count=count+1;
            end
        end
        
         y2(i+1)=y2(i)+0.5*h*(((l1*y2(i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i))))+((l1*y_2(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l1*(atan(x(i+1)))))));
         count=0;
         while(count<10)
            y2(i+1)=y2(i)+0.5*h*(((l1*y2(i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i))))+((l1*y2(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l1*(atan(x(i+1)))))));
            count=count+1;
         end
         
         y3(i+1)=y3(i)+0.5*h*(((l2*y3(i))+(1/(1+(x(i)*x(i))))-(l2*(atan(x(i))))+((l2*y_3(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l2*(atan(x(i+1)))))));
         count=0;
         while(count<10)
            y3(i+1)=y3(i)+0.5*h*(((l2*y3(i))+(1/(1+(x(i)*x(i))))-(l2*(atan(x(i))))+((l2*y3(i+1))+(1/(1+(x(i+1)*x(i+1))))-(l2*(atan(x(i+1)))))));
            count=count+1;
         end
         
        end
end
