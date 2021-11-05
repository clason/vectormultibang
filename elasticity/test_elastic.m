% TEST_ELASTIC test script for linearized elasticity examples
% This m-file computes discrete-valued optimal controls for the linearized
% elasticity equations using the approach and the parameters described in
%   C. Clason, C. Tameling, B. Wirth,
%   Vector-valued multibang control of differential equations
%   SICON 56 (2018), 2295-2326 / https://arxiv.org/abs/1611.07853
%
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

%% parameters
N = 65;             % number of vertices per dimension
E = 20;             % elastic modulus
nu = 0.3;           % Poisson's ratio
Om = [0,1,0,2];     % domain Omega: rectangle with corners (0,0) and (1,2)

% assemble finite element mesh, matrices
[Ah,Mh,xx,yy] = assembleElasticFEM(N,E,nu,Om);

%% choose examples
example = 4;
switch example
    case 1  % Figure 8a
        target  = 'rotation';
        penalty = 'radial';
        M       = 3;            % number of admissible states
        omega0  = 2*sqrt(2);    % magnitude of admissible states
        alpha   = 1e-3;         % multibang penalty parameter
    case 2  % Figure 8b
        target  = 'rotation';
        penalty = 'radial';
        M       = 5;            % number of admissible states
        omega0  = 2*sqrt(2);    % magnitude of admissible states
        alpha   = 1e-3;         % multibang penalty parameter
    case 3  % Figure 8c
        target  = 'rotation';
        penalty = 'concentric';
        alpha   = 1e-3;         % multibang penalty parameter
    case 4  % Figure 8d
        target  = 'rotation';
        penalty = 'concentric';
        alpha   = 1e-5;         % multibang penalty parameter
    case 5  % Figure 8e
        target  = 'attainable';
        penalty = 'concentric';
        alpha   = 1e-5;         % multibang penalty parameter
    case 6  % Figure 8d
        target  = 'perturbed';
        penalty = 'concentric';
        alpha   = 1e-5;         % multibang penalty parameter
end

%% construct target
switch target
    case 'rotation'     % rotate block phi degrees around center
        phi = -pi/6;
        center = [.5;1];
        B = [cos(phi) sin(phi);-sin(phi) cos(phi)]-eye(2);
        z = (B*[xx(:)-center(1) yy(:)-center(2)]')';
        z = z(:);
    case 'attainable'   % attainable target produced by force at beam end
        ue = zeros([size(xx),2]);
        ue(end,:,1) = 30;
        z = Ah\(Mh*ue(:));
    case 'perturbed'    % small random perturbation of attainable target
        ue = zeros([size(xx),2]);
        ue(end,:,1) = 30;
        z = Ah\(Mh*ue(:)) + .01*randn(2*N*N,1);
end

%% setup penalty
switch penalty
    case 'radial'
        ubangle = linspace(-pi,pi-2*pi/M,M);
        phi = ubangle + pi/M;     % calculate the midpoints between the angles
        phi(phi >= pi) = phi(phi >= pi) - 2*pi; % make sure that the angle is still between -pi and pi
        ubmag = omega0*ones(1,M);
        ubins = unique([-pi,phi,pi]);
        [ub1,ub2] = pol2cart(ubangle,ubmag);
        ub = [ub1;ub2];
        mb_penalty = @(p,gamma) mb_radial(p,alpha,gamma,ub,ubins,omega0);
    case 'concentric'
        ub = [1  1 -1 -1 2  2 -2 -2;...
              1 -1  1 -1 2 -2  2 -2];
        [ubangle,ubmag] = cart2pol(ub(1,:),ub(2,:));
        mb_penalty = @(p,gamma) mb_concentric(p,alpha,gamma,ub);
end
% color-coded plot of controls scaled to maximum admissible control value
uplot = @(u) phaseplot(xx,yy,u,ubmag);
% plot grid deformed according to state and target
yplot = @(y) deformplot(xx,yy,y,z);

%% solve control problem
[u,y] = ssn(z,Ah,Mh,N,mb_penalty,uplot,yplot);

%% visualize results
uplot(u);  yplot(y)
