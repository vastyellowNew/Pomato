function [jc] = jacobi(pos,vel,mup);
%This function calculates the Jacobi Constant for the circular restricted 3
%body problem state with pos=[x y z] and vel=[xd yd zd] in barycentered
%rotating frame coords for mass parameter mup.
    d=sqrt((pos(1)+mup)^2+pos(2)^2+pos(3)^2);
    r=sqrt((pos(1)-1.0+mup)^2+pos(2)^2+pos(3)^2);
    ustar=(1.0-mup)/d+mup/r+.5*(pos(1)^2+pos(2)^2);
    vsqr=vel(1)^2+vel(2)^2+vel(3)^2;
    jc=2*ustar-vsqr;
return
