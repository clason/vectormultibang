function [A,M,xx,yy] = assembleElasticFEM(N,E,nu,Om)
% ASSEMBLEELASTICFEM assemble finite element mesh and matrices for linearized elasticity
% [A,M,XX,YY] = ASSEMBLEELASTICFEM(N,E,NU,OM) assembles the structured 
% finite element mesh (returned as vector of vertex coordinates [XX,YY]), 
% stiffness matrix A, and mass matrix M for the equations of linearized 
% elasticity. N is the number of vertices per dimension, E and NU are the
% material constants (elastic modulus and Poisson ratio, respectively), and
% OM = [ax ay bx by] specifies the dimensions of the computational domain,
% which is given by a box with opposite vertices (ax,ay) and (bx,by).
% Dirichlet conditions on the bottom boundary xx == ax are prescribed.
%     
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

%% setup mesh
[xx,yy] = meshgrid(linspace(Om(1),Om(2),N),linspace(Om(3),Om(4),N));
nodes   = [xx(:) yy(:)];
numNodes = length(nodes);

%% assemble stiffness, mass matrix
elements = delaunay(nodes(:,1),nodes(:,2));
numElems = size(elements,1);

rows = zeros(36*length(elements),1);
cols = rows;
valsM = rows;
valsL = rows;
valsK = rows;

I = eye(2);
D = [-1 -1;I]';

counter = 1;
for k = 1:numElems                    % loop through all elements {T}
    nodeInds = elements(k, 1:3);      % get global indices
    x = nodes(nodeInds, 1:2);         % global nodes (as matrix rows)
    A = (D*x)';                       % Jacobian of coordinate transform from T_ref to T.
    AinvD = A'\D;                     % gradient operator
    div = [[1 0]*AinvD [0 1]*AinvD];  % divergence operator
    eps = [AinvD zeros(2,3);          % symmetrized gradient operator
           zeros(2,3) AinvD];
    eps(2,:) = (eps(2,:)+eps(3,:))/2;
    eps(3,:) = eps(2,:);              
    vol = abs(det(A)/2);              % element volume
    
    % assemble matrices
    [c,r] = meshgrid([nodeInds nodeInds+numNodes],[nodeInds nodeInds+numNodes]);
    rows(counter:counter+35) = r(:);
    cols(counter:counter+35) = c(:);
    valsM(counter:counter+35) = [[2 1 1;1 2 1;1 1 2] zeros(3);
                                  zeros(3) [2 1 1;1 2 1;1 1 2]]*(vol/12);
    valsL(counter:counter+35) = vol*(eps'*eps);
    valsK(counter:counter+35) = vol*(div'*div);
    counter = counter + 36;
end

M = sparse(rows,cols,valsM,2*numNodes,2*numNodes);
L = sparse(rows,cols,valsL,2*numNodes,2*numNodes);
K = sparse(rows,cols,valsK,2*numNodes,2*numNodes);

%% differential operator
A = 2*E/(2*(1+nu))*L + E*nu/((1+nu)*(1-2*nu))*K;

%% implement Dirichlet data
dirichletNodes = 0*xx;
dirichletNodes(1,:) = 1;
dirichletNodes = logical(dirichletNodes);
M([dirichletNodes(:);dirichletNodes(:)],:) = 0;
A([dirichletNodes(:);dirichletNodes(:)],:) = 0;
A(:,[dirichletNodes(:);dirichletNodes(:)]) = 0;
for i = 1:numNodes
    if dirichletNodes(i)
        A(i,i) = 1;
        A(i+numNodes,i+numNodes) = 1;
    end
end