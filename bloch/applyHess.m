function Hdx = applyHess(d,Xk,dx)
% APPLYHESS compute application of Hessian
% HDX = APPLYHESS(D,XK,DX) computes the action HDX of the Hessian in 
% direction DX. The structure XK contains the point in which the Hessian
% is evaluated. The structure D contains the problem parameters.
%     
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)
% adapted from https://github.com/chaigner/rfcontrol
% February 3, 2015         Christoph S. Aigner (christoph.aigner@tugraz.at)
%                          Christian Clason (christian.clason@uni-due.de)
%                          Armin Rund (armin.rund@uni-graz.at)

N = Xk.N;     P = Xk.P;    u = Xk.u;  v = Xk.v;    w = d.w;

du = dx(1:d.Nu); dv = dx(1+d.Nu:end);

steps = length(d.tdis);
dt    = d.dt;
Nu    = d.Nu;

I     = eye(3,3);
B1    = d.gyro*d.B1c;
B3    = d.gyro*d.G3;

%% compute directional derivatives of state, adjoint
dN = zeros(3,d.Nw,steps-1);
dP = zeros(3,d.Nw,steps-1);
parfor z = 1:d.Nw
  Nz = N(:,z,:); Pz = P(:,z,:);

  % solve linearized state equation
  dMz  = zeros(3,steps);
  for k = 2:steps
      bk = [ B1*Nz(3,k-1)*dv(k-1);...
             B1*Nz(3,k-1)*du(k-1); ...
            -B1*Nz(2,k-1)*du(k-1) - B1*Nz(1,k-1)*dv(k-1)];
      Ak = [               0,     w(k-1,z)*B3,       v(k-1)*B1;...
                -w(k-1,z)*B3,               0,       u(k-1)*B1;...
                  -v(k-1)*B1,      -u(k-1)*B1,               0 ];
      dMz(:,k) = (I-dt/2*Ak) \ ((I+dt/2*Ak)*dMz(:,k-1) + dt*bk);
  end
  dN(:,z,:) = (dMz(:,1:Nu)+dMz(:,2:Nu+1))/2;    %projection of state onto control space
  
  % solve linearized adjoint equation
  Akp1 = Ak';
  bkp1 = [-B1*Pz(3,end)*dv(end);...
          -B1*Pz(3,end)*du(end);... 
           B1*Pz(2,end)*du(end) + B1*Pz(1,end)*dv(end)];
  dPz  = zeros(3,steps-1);
  dPz(:,end) = (I-dt/2*Akp1)\(-dMz(:,end)+dt/2*bkp1);
  for k = steps-2:-1:1
      bk = [-B1*Pz(3,k)*dv(k);...
            -B1*Pz(3,k)*du(k);... 
             B1*Pz(2,k)*du(k) + B1*Pz(1,k)*dv(k)];
      Ak = [               0,      -w(k,z)*B3,        -v(k)*B1;...
                   w(k,z)*B3,               0,        -u(k)*B1;...
                     v(k)*B1,         u(k)*B1,               0 ];
      dPz(:,k) = (I-dt/2*Ak) \ ((I+dt/2*Akp1)*dPz(:,k+1) + dt/2*(bk+bkp1));
      Akp1 = Ak;    bkp1 = bk;
  end
  dP(:,z,:) = dPz;
    
end;

% action of Hessian (directional derivatives of state)
dQ1 = B1*squeeze(sum(dN(3,:,:).* P(2,:,:) - dN(2,:,:).* P(3,:,:) + ...
                      N(3,:,:).*dP(2,:,:) -  N(2,:,:).*dP(3,:,:),2));
dQ2 = B1*squeeze(sum(dN(3,:,:).* P(1,:,:) - dN(1,:,:).* P(3,:,:) + ...
                      N(3,:,:).*dP(1,:,:) -  N(1,:,:).*dP(3,:,:),2));

Hdx = dx - Xk.DHg*[dQ1;dQ2];
