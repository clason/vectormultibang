% TEST_BLOCH test script for Bloch examples
% This m-file computes discrete-valued optimal controls for the Bloch
% equation using the approach and the parameters described in
%   C. Clason, C. Tameling, B. Wirth,
%   Vector-valued multibang control of differential equations
%   https://arxiv.org/abs/161X.XXXXXX
%     
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

%% set default parameters
% time discretization
d.T     = 7;                      % optimization time in ms
d.Nt    = 1001;                   % total number of temporal points
d.tdis  = linspace(0,d.T,d.Nt);   % temporal running variable
d.dt    = d.tdis(2) - d.tdis(1);  % temporal grid size
d.Nu    = d.Nt-1;                 % number of temporal control points

% model parameters
d.gyro  = 267.51;                 % gyromagnetic ratio 42.57*2*pi in MHz/T
d.B0    = 3000;                   % static magnetic field strength in mT
d.M0c   = 1;                      % normalized equilibrium magnetization
d.B1c   = 1e-2;                   % weighting for the RF amplitude [u*1e3*d.B1c] = muT
d.G3    = 1;                      % weighting for the z-Gradient in mT
d.u0 = zeros(d.Nu,1);             % RFx initial guess
d.v0 = zeros(d.Nu,1);             % RFy initial guess
d.w  = ones(d.Nu,1);              % fixed with external shape

% problem parameters
d.alpha  = 1e-1;                  % control costs for u (SAR)
d.M0     = [0;0;1];               % initial magnetization
d.Md     = [1;0;0];               % target  magnetization
d.M      = 6;                     % number of control states
d.omega0 = 1;                     % radius
d.theta0 = 0;                     % offset (must be smaller than 2*pi/d.M)

%% vary parameters for different examples
example = 4;
switch example
    case 1 % Figure 4
        d.w = d.w*1e-2;
        d.M = 3;
    case 2 % Figure 5
        d.w = d.w*1e-2;
        d.M = 6;
    case 3 % Figure 6
        d.w = d.w*1e-2*[1,2,3,4];
        d.M = 6;
        d.Md = [d.Md d.Md d.Md d.Md];
    case 4 % Figure 7
        d.w = d.w*1e-2*[1,2,3,4];
        d.M = 6;
        d.Md = [d.M0 d.M0 d.Md d.M0];
end

% compute admissible control states
d.Nw = size(d.w,2);
d.ub = zeros(2,d.M);          % initialize vector for control states
d.ubangle = d.theta0 + linspace(-pi,pi-2*pi/d.M,d.M);
d.phi = d.ubangle + pi/d.M;   % calculate the midpoints between the angles
d.phi(d.phi >= pi) = d.phi(d.phi >= pi) - 2*pi; % make sure that the angle is still between -pi and pi
d.ubmag = d.omega0*ones(1,d.M);
d.ubins = unique([-pi,d.phi,pi]);
[ub1,ub2] = pol2cart(d.ubangle,d.ubmag);
d.ub = [ub1;ub2];

%% solve control problem
x = ssn(d,@objfun,@applyHess);
    
%% plot solution
blochplot(d,x)