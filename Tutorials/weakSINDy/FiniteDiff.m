%% Finite difference
function Dev=FiniteDiff(x,dx,d)

[n,m]=size(x);

[len,Index]=max([n,m]);

if Index==1
    Dev=zeros(len,1);
else
    Dev=zeros(1,len);
end

if d==1
    for i=2:len-1
        Dev(i)=(x(i+1)-x(i-1))/(2*dx);
    end
    Dev(1) = (-3/2*x(1) + 2*x(2) - x(3)/2)/dx;
    Dev(end) = (3.0/2*x(end) - 2*x(end-1) + x(end-2)/2)/dx;
elseif d==2
    for i=2:len-1
        Dev(i)=(x(i+1)-2*x(i)+x(i-1)) / dx^2;
    end
    Dev(1) = (2*x(1) - 5*x(2) + 4*x(3) - x(4)) / dx^2;
    Dev(end) = (2*x(end) - 5*x(end-1) + 4*x(end-2) - x(end-3)) / dx^2;
elseif d==3
    for i=3:len-2
        Dev(i) = (x(i+2)/2-x(i+1)+x(i-1)-x(i-2)/2) / dx^3;
    end
        
    Dev(1) = (-2.5*x(1)+9*x(2)-12*x(3)+7*x(4)-1.5*x(5)) / dx^3;
    Dev(2) = (-2.5*x(2)+9*x(3)-12*x(4)+7*x(5)-1.5*x(6)) / dx^3;
    
    Dev(end) = (2.5*x(end)-9*x(end-1)+12*x(end-2)-7*x(end-3)+1.5*x(end-4)) / dx^3;
    Dev(end-1) = (2.5*x(end-1)-9*x(end-2)+12*x(end-3)-7*x(end-4)+1.5*x(end-5)) / dx^3;
elseif d>3
    Dev=FiniteDiff(FiniteDiff(x,dx,3), dx, d-3);
end







