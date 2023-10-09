function dx = roessler(t,x,p)

dx = [
    -x(2)-x(3); 
    x(1)+p(1)*x(2); 
    p(2)+(x(1)-p(3))*x(3)
    ];