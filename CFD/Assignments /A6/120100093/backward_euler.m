function [y1,y2,y3] = backward_euler(a,y_1,y_2,y_3,x,h,n,N,l1,l2)

y1(1)=a;
y2(1)=a;
y3(1)=a;
err=1;
epsilon=1*10^-6;

for i=2:n+1
    if N==0
        y1(i)=y1(i-1)+h*((x(i)*cos(x(i)))+(y_1(i)/x(i)));
        err=y1(i)-y_1(i);
        count=0;
        while(abs(err)>epsilon && count<100)
            temp=y1(i);
            y1(i)=y1(i-1)+h*((x(i)*cos(x(i)))+(y1(i)/x(i)));
            err=y1(i)-temp;
            count=count+1;
        end
        
        y2(i)=y2(i-1)+h*((l1*y_2(1,i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i)))));
        err=y2(i)-y_2(i);
        count=0;
        while(abs(err)>epsilon && count<100)
            temp=y2(i);
            y2(i)=y2(i-1)+h*((l1*y2(1,i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i)))));
            err=y2(i)-temp;
            count=count+1;
        end
        
        y3(i)=y3(i-1)+h*((l2*y_3(1,i))+(1/(1+(x(i)*x(i))))-(l2*atan(x(i))));
        err=y3(i)-y_3(i);
        count=0;
        while(abs(err)>epsilon && count<100)
            temp=y3(i);
            y3(i)=y3(i-1)+h*((l2*y3(1,i))+(1/(1+(x(i)*x(i))))-(l2*(atan(x(i)))));
            err=y3(i)-temp;
            count=count+1;
        end
    end
    
    if N==1
         y1(i)=y1(i-1)+h*((x(i)*cos(x(i)))+(y_1(i)/x(i)));
         y1(i)=y1(i-1)+h*((x(i)*cos(x(i)))+(y1(i)/x(i)));
         
         y2(i)=y2(i-1)+h*((l1*y_2(1,i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i)))));
         y2(i)=y2(i-1)+h*((l1*y2(1,i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i)))));
         
         y3(i)=y3(i-1)+h*((l2*y_3(1,i))+(1/(1+(x(i)*x(i))))-(l2*atan(x(i))));
         y3(i)=y3(i-1)+h*((l2*y3(1,i))+(1/(1+(x(i)*x(i))))-(l2*(atan(x(i)))));
    end
    
    if N==10
         y1(i)=y1(i-1)+h*((x(i)*cos(x(i)))+(y_1(i)/x(i)));
         count=0;
         while(count<10) 
            y1(i)=y1(i-1)+h*((x(i)*cos(x(i)))+(y1(i)/x(i)));
            count=count+1;
         end
            
         y2(i)=y2(i-1)+h*((l1*y_2(1,i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i)))));
         count=0;
         while(count<10)
            y2(i)=y2(i-1)+h*((l1*y2(1,i))+(1/(1+(x(i)*x(i))))-(l1*(atan(x(i)))));
            count=count+1;
         end
         
         y3(i)=y3(i-1)+h*((l2*y_3(1,i))+(1/(1+(x(i)*x(i))))-(l2*atan(x(i))));
         count=0;
         while(count<10)
            y3(i)=y3(i-1)+h*((l2*y3(1,i))+(1/(1+(x(i)*x(i))))-(l2*(atan(x(i)))));
            count=count+1;
         end
    end
end
