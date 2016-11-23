function P = cn_adjoint(d,PT,u,v,w)
% CN_ADJOINT solve adjoint Bloch equation using adjoint Crank-Nicolson scheme
% M = CN_ADJOINT(D,PT,U,V,W) computes the magnetization vector M starting from
% terminal conditions PT with RF pulse U,V and gradient W. The structure D
% contains the problem parameters.
%
% adapted from https://github.com/chaigner/rfcontrol
% February 3, 2015         Christoph S. Aigner (christoph.aigner@tugraz.at)
%                          Christian Clason (christian.clason@uni-due.de)
%                          Armin Rund (armin.rund@uni-graz.at)

steps = length(d.tdis);
dt    = d.dt;

I  = eye(3,3);
B1 = d.gyro*d.B1c;
B3 = d.gyro*d.G3;

P = zeros(3,d.Nw,steps-1);
parfor z = 1:d.Nw
    Akp1 = [               0,    -w(end,z)*B3,      -v(end)*B1;...
                 w(end,z)*B3,               0,      -u(end)*B1;...
                   v(end)*B1,       u(end)*B1,               0 ];
    Pz = zeros(3,steps-1);
    Pz(:,end) = (I-dt/2*Akp1)\PT(:,z,:);
    for k = steps-2:-1:1
        Ak = [               0,       -w(k,z)*B3,        -v(k)*B1;...
                     w(k,z)*B3,                0,        -u(k)*B1;...
                       v(k)*B1,          u(k)*B1,               0 ];
        Pz(:,k) = (I-dt/2*Ak)\((I+dt/2*Akp1)*Pz(:,k+1));
        Akp1 = Ak;
    end
    P(:,z,:) = Pz;
end