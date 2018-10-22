%Fixed-Point Computation Tests
%A J-E Planar Map Test Case
%Author - Wayne Schlei
%Date - 11/29/12
clear all; close all; clc;

%Numerical Simulation
tol=1e-12;
options=odeset('RelTol',tol,'AbsTol',tol,'Events',@positiveycross);
options2=odeset('RelTol',tol,'AbsTol',tol);
options3=odeset('RelTol',tol,'AbsTol',tol,'Events',@stm2Devent);
options4=odeset('RelTol',tol,'AbsTol',tol,'Events',@posycrossnoterm);
optionsL5 = odeset('RelTol',tol,'AbsTol',tol,'Events',@xL5);

%% System Parameters
mup = 3202.739 /(3202.739 + 1.26686535e8)

%Pick a Jacobi constant
jc=3.00
%Pick Number of Crossings
numCross = 1
%xGuess = [-1.205, 0.0]; numCross = 1;
%xGuess = [-1.268, 0.0]; numCross = 3;
xGuess = [-1.123, 0.0]; numCross = 1;
xGuess = [-1.154, 0.0]; numCross = 1;
xGuess = [-1.305, 0.0]; numCross = 1;
xGuess = [-1.24, 0.0]; numCross = 1;

%IC
y= 0.0;
x= xGuess(1); xd= xGuess(2);
%True fixed point at (-1.20362,0)
d=sqrt((x+mup)^2+y^2);
r=sqrt((x-1+mup)^2+y^2);
ustar=(1-mup)/d+mup/r+.5*(x^2+y^2);
yd=sqrt(2*ustar-jc-xd^2);
stm0 = eye(6,6);

%Plot of whole trajectory
figure(1);
plot(-mup,0,'ro');
hold on;
plot(1-mup,0,'o','color',[.4 .4 .4]);
axis([-1.5 1.5 -1.5 1.5]);
axis equal;
xlabel('\itx\rm (Nondim)');
ylabel('\ity\rm (Nondim)');
title('Jupiter-Europa Test Case -WS');

% A simple test of STM
% x0=[x,y,0,xd,yd,0];
% nom = x0;    
% x0=[nom,stm0(1,:),stm0(2,:),stm0(3,:),stm0(4,:),stm0(5,:),stm0(6,:)];
% [tau,xsol] = ode45(@(tau,xsol) stm(tau,xsol,mup),[0 2.05],x0,options2);
% %Calculate STM
% phi=reshape(xsol(end,7:end),6,6)'

%% Map Computation
if (isreal(yd)==1)
    x0=[x,y,0,xd,yd,0];
    nom = x0;
    for k=1:numCross
        x0=[nom,stm0(1,:),stm0(2,:),stm0(3,:),stm0(4,:),stm0(5,:),stm0(6,:)];
        dtau=[0 1000];
        %[tau,xsol]=ode45(@(tau,xsol) stm(tau,xsol,mup),dtau,x0,options);
        [tau,xsol]=ode45(@(tau,xsol) cr3bp(tau,xsol,mup),dtau,nom,options);
        %Save the x and xd points of the crossings
        x_save(k)=xsol(end,1);
        xd_save(k)=xsol(end,4);
        %Grab the previous crossing as the next crossing
        nom=xsol(end,1:6);
        %Plot the crossings
        plot(nom(1),nom(2),'g.');
        plot(xsol(:,1),xsol(:,2),'b');
    end
end
period = tau(end)

% %% fsolve to find fixed point
% x0 = [-1.205; 0];           % Make a starting guess at the solution
% %x0 = [-1.20362; 0]; %The answer...
% opts=optimset('Display','iter','Jacobian','on');   % Option to display output
% [xsolve,fval] = fsolve(@(x) fixedPointFunc(x,jc,mup,numCross),x0,opts)  % Call optimizer

% %% STD Newton
% eps = 1;
% x = x0;
% dx = [0;0];
% while (eps>1e-10)
%     x = x + dx
%     [F, J] = fixedPointFunc(x,jc,mup,1);
%     eps = norm(F)
%     dx = - J'*inv(J*J')*F;
% end

%% Multiple-shooting with Periodicity, Jacobi Constant, and Section Consts
k = 8*numCross
T = period;
%Subsample trajectory
dt = [0:(k-1)]*T/k;
state = [x,y,0,xd,yd,0];
[tau,xsol]=ode45(@(tau,xsol) cr3bp(tau,xsol,mup),dt,state,options2);
X = [];
for i=1:k
    X = [X;xsol(i,:)'];
end
X = [X;T]; %Free-variable vector with 6*k+1 elements

Fmag = 1000;
eps = 1.e-12;
dX = zeros(6*k+1,1);
F = zeros(6*k, 1);
DF = zeros(6*k, 6*k+1);
for i=1:k-1
	DF( 6*(i-1)+1:6*i , 6*i+1:6*(i+1) ) = -eye(6,6);
end
%Psi_1 - Periodicity Constraint derivatives wrt x0
DF(6*(k-1)+1, 1) = -1;
DF(6*(k-1)+2, 2) = -1;
DF(6*(k-1)+3, 3) = -1;
DF(6*(k-1)+4, 4) = -1;
DF(6*(k-1)+5, 6) = -1;
DF(6*k,2) = 1; %Section constraint -> y=0

%Newton Solution loop
maxIter = 20; iter = 0;
while (Fmag >= eps)
	if (iter >= maxIter)
		break;
	end
	
	%Integrate Arcs
	for i=1:k
		nom = X(6*(i-1)+1:6*(i-1)+6)';
		IC = [nom,stm0(1,:),stm0(2,:),stm0(3,:),stm0(4,:),stm0(5,:),stm0(6,:)];
		dt = [0 X(6*k+1)/k];
		[tau,xsol]=ode45(@(tau,xsol) stm(tau,xsol,mup),dt,IC,options2);
		%Final State information (t_f)
		xf=xsol(end,1);yf=xsol(end,2);zf=xsol(end,3);
		xdf=xsol(end,4);ydf=xsol(end,5);zdf=xsol(end,6);
		dstate=cr3bp(tau(end),[xf,yf,zf,xdf,ydf,zdf],mup);
		%Calculate STM
		phi=reshape(xsol(end,7:end),6,6)';
		%Compute F(X)
		if (i<=(k-1))
			F(6*(i-1)+1:6*(i-1)+6) = (xsol(end,1:6) - X(6*i+1:6*i+6)')';
			%Input STM's into DF(X)
			DF(6*(i-1)+1:6*(i-1)+6,6*(i-1)+1:6*(i-1)+6) = phi;
			%Inpute state derivatives into DF(X)
			DF(6*(i-1)+1:6*(i-1)+6,6*k+1) = dstate'/k;
		else
			%i=k case -> Periodicity
            %5 Equations in F
			F(6*(i-1)+1:6*(i-1)+4) = (xsol(end,1:4) - X(1:4)')';
			%skip ydot constraint
			F(6*(i-1)+4) = xsol(end,6)-X(6);
            
			%STM elements 5x6
			DF(6*(i-1)+1:6*(i-1)+4,6*(i-1)+1:6*(i-1)+6) = phi(1:4,:);
			DF(6*(i-1)+4, 6*(i-1)+1:6*(i-1)+6) = phi(6,:);%Skip ydot
			%Input State Derivatives
            DF(6*(i-1)+1:6*(i-1)+4,6*k+1) = dstate(1:4)'/k;
			DF(6*(i-1)+4,6*k+1) = dstate(6)/k;%Skip ydot
		end
		if (i == 1)
            %Section Constraint equation
            F(6*k) = X(2);
            DF(6*k,1:6) = [0 1 0 0 0 0];
			%Jacobi constant constraint
			F(6*k+1) = jacobi(nom(1:3),nom(4:6),mup) - jc;
			%DF matrix for Jacobi
			Ui = cr3bpPartials(nom,mup);
			DF(6*k+1,1:6) = 2*[Ui,-1*nom(4:6)];  
        end
    end
    %Update info
    sprintf('iter = %d',iter)
    Fmag = norm(F);
    sprintf(' norm(F) = %6.4E ', Fmag)
    %Compute Update
    dX = -inv(DF)*F; %Newton (m=n)
    %dX = -DF'*inv(DF*DF')*F; %Min-Norm (m<n)
    X = X + dX;
    iter = iter + 1;
%     keyboard
end
if (iter >= maxIter)
	sprintf('Corrections terminated with too many iterations')
end

%Draw corrected orbit on figure(1)
for i=1:k
    nom = X(6*(i-1)+1:6*(i-1)+6);
    dtau=[0 X(6*k+1)/k];
    [tau,xsol]=ode45(@(tau,xsol) cr3bp(tau,xsol,mup),dtau,nom,options2);
    %Plot the segments
    plot(nom(1),nom(2),'r.');
    plot(xsol(:,1),xsol(:,2),'r');
end