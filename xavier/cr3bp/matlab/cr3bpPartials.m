function [first,second] = cr3bpPartials(state,mup);
%Evaluates the partial derivatives of Ustar given state and mup
x = state(1);
y = state(2);
z = state(3);
xd = state(4);
yd = state(5);
zd = state(6);

d=sqrt((x+mup)^2+y^2+z^2);
r=sqrt((x-1+mup)^2+y^2+z^2);

d3 = d*d*d; d5 = d*d*d3;
r3 = r*r*r; r5 = r*r*r3;

%First order partials
Ux = x - (1.0-mup)*(x+mup)/d3 - mup*(x-1.0+mup)/r3;
Uy = y - (1.0-mup)*y/d3 - mup*y/r3;
Uz = -(1.0-mup)*z/d3 - mup*z/r3;
first = [Ux,Uy,Uz];

if nargout > 1
	%Return the second order partials
	Uxx =1-(1-mup)/d3-mup/r3+3*(1-mup)*(x+mup)^2/d5+3*mup*(x-1+mup)^2/r5;
	Uxy =3*(1-mup)*(x+mup)*y/d5+3*mup*(x-1+mup)*y/r5;
	Uxz =3*(1-mup)*(x+mup)*z/d5+3*mup*(x-1+mup)*z/r5;
	Uyy =1-(1-mup)/d3-mup/r3+3*(1-mup)*y^2/d5+3*mup*y^2/r5;
	Uyz =3*(1-mup)*y*z/d5+3*mup*y*z/r5;
	Uzz =-(1-mup)/d3-mup/r3+3*(1-mup)*z^2/d5+3*mup*z^2/r5;
	second = [Uxx,Uxy,Uxz,Uyy,Uyz,Uzz];
end