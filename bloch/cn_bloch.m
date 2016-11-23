function M = cn_bloch(d,M0,u,v,w)
% CN_BLOCH solve Bloch equation using Crank-Nicolson scheme
% M = CN_BLOCH(D,M0,U,V,W) computes the magnetization vector M starting
% from initial conditions M0 with RF pulse U,V and gradient W. The
% structure D contains the problem parameters.
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

M = zeros(3,d.Nw,steps); 

parfor z = 1:d.Nw
    Mz = zeros(3,steps); Mz(:,1) = M0;
    for k = 2:steps
        Ak = [               0,       w(k-1,z)*B3,       v(k-1)*B1;...
                  -w(k-1,z)*B3,                 0,       u(k-1)*B1;...
                    -v(k-1)*B1,        -u(k-1)*B1,               0 ];
        Mz(:,k) = (I-dt/2*Ak)\((I+dt/2*Ak)*Mz(:,k-1));
    end
    M(:,z,:) = Mz;
end