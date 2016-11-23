function [J,G,Xk] = objfun(d,x)
% OBJFUN compute functional value, gradient
% [J,G,XK] = OBJFUN(D,X) computes the value J of the functional 
% to be minimized together with the gradient G in the point X=(u,v).
% The structure XK contains the necessary information to evaluate the 
% Hessian in X. The structure D contains the problem parameters.
%     
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

u = x(1:d.Nu);  v = x(1+d.Nu:end);

%% objective value
M = cn_bloch(d,d.M0,u,v,d.w);   % state 
res = d.Md - M(:,:,end);        % residual
ang = angle(u+1j*v);
bins = unique([-pi,d.ubangle,pi]);
[~,~,i1] = histcounts(ang+1e-8,bins);   i1(i1==d.M+1)=1; i1(i1==0)=d.M; % wrap around
J = norm(res(:))^2/2 + d.gamma/2*d.dt*(norm(u)^2+norm(v)^2) ...
    + d.alpha/2*d.ubmag(1)*d.dt*(cos(d.phi(i1))*u+sin(d.phi(i1))*v)/sin((d.ubangle(2)-d.ubangle(1))/2);

%% gradient
P = cn_adjoint(d,res,u,v,d.w);
N = 1/2*(M(:,:,1:end-1) + M(:,:,2:end)); % projection of state onto control space (from piecewise linear to piecewise constant)
B1 = d.gyro*d.B1c;
Q1 = B1*squeeze(sum(N(3,:,:).*P(2,:,:) - N(2,:,:).*P(3,:,:),2));
Q2 = B1*squeeze(sum(N(3,:,:).*P(1,:,:) - N(1,:,:).*P(3,:,:),2));

% compute active sets, regularized subdifferential, Newton derivative
[Hg,DHg,as,rnodes] = mb_radial([Q1,Q2],d.alpha,d.gamma,d.ub,d.ubins,d.omega0);

G = [u;v] - [Hg(:,1); Hg(:,2)];

%% Hessian information
Xk.N = N;
Xk.P = P;
Xk.u = u;
Xk.v = v;
Xk.Q1 = Q1;
Xk.Q2 = Q2;
Xk.DHg = DHg;
Xk.As = as;
Xk.rnodes = rnodes;