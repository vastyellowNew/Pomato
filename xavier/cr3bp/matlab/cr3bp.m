function xdot=cr3bp(tau,xh,mup)
%Author - Wayne Schlei    Mod. Date - 1/21/08
%  Equations of Motion corresponding to the CR3BP in the
%rotating frame and all quantities are nondimensional 
%such that nondim G=1 and nondim n=1.  Note: derivatives are with
%respect to the nondimensional time parameter tau.
% %Mass Parameter - Earth-Moon
% %Earth-Moon System mu=.01215
% mum=4902.7949;mue=3.986004418e5;
% mup=mum/(mum+mue);
%Setting Variables
x=xh(1);
y=xh(2);
z=xh(3);
xd=xh(4);
yd=xh(5);
zd=xh(6);
%Relative Velocities (First 3 first order differential equations)
xdot(1)=xd;
xdot(2)=yd;
xdot(3)=zd;
%Relative distances
d=sqrt((x+mup)^2+y^2+z^2);
r=sqrt((x-1+mup)^2+y^2+z^2);
%Last 3 first order differential equations (From EOM's)
xdot(4)=2*yd+x-(1-mup)*(x+mup)/(d^3)-mup*(x-1+mup)/(r^3);
xdot(5)=-2*xd+y-(1-mup)*y/(d^3)-mup*y/(r^3);
xdot(6)=-(1-mup)*z/(d^3)-mup*z/(r^3);
xdot=xdot';
return
