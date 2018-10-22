function xdot=stm(tau,xh,mup)
%Author - Wayne Schlei    Mod. Date - 2/26/08
%  This program is the first order differential equations to find the STM
%for the CR3BP in full 3D space (x,y,z).  xh is a 42 element vector with 
%the reference solution states as the first 6 entries and the elements of
%the STM as the final 36 entries.
%=========================================================================
% %Mass Parameter - Earth-Moon
% %Earth-Moon System mu=.01215
% mum=4902.7949;mue=3.986004418e5;
% mup=mum/(mum+mue);
%=================================================================
%Reference Solution
%  6 scalar Equations to evalutate the reference solution
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
%==================================================================
%A(t) Matrix
A11=zeros(3,3);
A12=eye(3,3);
%Partials evaluated at reference solution
dfxdx=1-(1-mup)/(d^3)-mup/(r^3)+3*(1-mup)*(x+mup)^2/(d^5)+3*mup*(x-1+mup)^2/(r^5);
dfxdy=3*(1-mup)*(x+mup)*y/(d^5)+3*mup*(x-1+mup)*y/(r^5);
dfydx=dfxdy;
dfxdz=3*(1-mup)*(x+mup)*z/(d^5)+3*mup*(x-1+mup)*z/(r^5);
dfzdx=dfxdz;
dfydy=1-(1-mup)/(d^3)-mup/(r^3)+3*(1-mup)*y^2/(d^5)+3*mup*y^2/(r^5);
dfydz=3*(1-mup)*y*z/(d^5)+3*mup*y*z/(r^5);
dfzdy=dfydz;
dfzdz=-(1-mup)/(d^3)-mup/(r^3)+3*(1-mup)*z^2/(d^5)+3*mup*z^2/(r^5);
A21=[dfxdx,dfxdy,dfxdz;dfydx,dfydy,dfydz;dfzdx,dfzdy,dfzdz];
A22=[0,2,0;-2,0,0;0,0,0];
A=[A11,A12;A21,A22];
%==================================================================
%STM phi(t,t0)
%    36 elements of the matrix to give solution of variational equations
%    dxf=phi*dx0
%Define Elements
% phi(1,1)=xh(7);phi(1,2)=xh(8);phi(1,3)=xh(9);phi(1,4)=xh(10);phi(1,5)=xh(11);phi(1,6)=xh(12);
% phi(2,1)=xh(13);phi(2,2)=xh(14);phi(2,3)=xh(15);phi(2,4)=xh(16);phi(2,5)=xh(17);phi(2,6)=xh(18);
% phi(3,1)=xh(19);phi(3,2)=xh(20);phi(3,3)=xh(21);phi(3,4)=xh(22);phi(3,5)=xh(23);phi(3,6)=xh(24);
% phi(4,1)=xh(25);phi(4,2)=xh(26);phi(4,3)=xh(27);phi(4,4)=xh(28);phi(4,5)=xh(29);phi(4,6)=xh(30);
% phi(5,1)=xh(31);phi(5,2)=xh(32);phi(5,3)=xh(33);phi(5,4)=xh(34);phi(5,5)=xh(35);phi(5,6)=xh(36);
% phi(6,1)=xh(37);phi(6,2)=xh(38);phi(6,3)=xh(39);phi(6,4)=xh(40);phi(6,5)=xh(41);phi(6,6)=xh(42);
counter=6;
for n=1:1:6;
    for k=1:1:6;
        counter=counter+1;
        phi(n,k)=xh(counter);
    end;
end;
%Differential Equations
phidot=A*phi;
%Populate Output matrix
counter=6;
for n=1:1:6;
    for k=1:1:6;
        counter=counter+1;
        xdot(counter)=phidot(n,k);
    end;
end;
xdot=xdot';
return