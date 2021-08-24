% TEST_BRANCHEDTRANSPORT test script for multimaterial branched transport
% This m-file computes discrete-valued optimal transport flows for the
% multimaterial branched transport using the approach and the parameters
% described in
%   C. Clason, C. Tameling, B. Wirth,
%   Convex relaxation of discrete vector-valued optimization problems
%   https://arxiv.org/abs/2108.10077
%   To appear in SIAM Review (2021)
%
% April 19, 2021                    Christian Clason (c.clason@uni-graz.at)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

%% generate transportation problem

example = 4;

% generate random street network
rng(12345);
x = linspace(0,1,10);
[X,Y] = meshgrid(x,x);
X = [.8;.5;.2;X(:);.2;.5;.8];
Y = [.1;0;.2;Y(:);.8;1;.9];
X = X + .03*randn(size(X));
Y = Y + .03*randn(size(Y));
param.N = length(X);
param.M = 3;  % number of transported materials
switch example
    case 1  % default - change nothing
    case 2  % swap third source and sink
        aux = X(1); X(1) = X(end-2); X(end-2) = aux;
        aux = Y(1); Y(1) = Y(end-2); Y(end-2) = aux;
    case 3  % swap first and third sink
        aux = X(end); X(end) = X(end-2); X(end-2) = aux;
        aux = Y(end); Y(end) = Y(end-2); Y(end-2) = aux;
    case 4  % add fourth material
        param.M = 4;
end
param.nodes = [X,Y];
DT = delaunayTriangulation(param.nodes);
triplot(DT);
param.edges = edges(DT);
param.lengths = sqrt(sum((param.nodes(param.edges(:,1),:)-param.nodes(param.edges(:,2),:)).^2,2));
param.edges(param.lengths>.2,:) = [];
param.numE = length(param.edges);
param.lengths = sqrt(sum((param.nodes(param.edges(:,1),:)-param.nodes(param.edges(:,2),:)).^2,2));

% generate sources and sinks
masses = ones(param.M,1);
sources = zeros(param.N,param.M);
sinks = sources;
for j = 1:param.M
    sources(j,j) = masses(j);
    sinks(end+1-j,end+1-j) = masses(j);
end
param.source = sources - sinks;
hold on;
plot(X([1:param.M,end+1-param.M:end]),Y([1:param.M,end+1-param.M:end]),'ko');
axis equal;
hold off;

% generate multibang penalty
ub = zeros(param.M,2^param.M);
for j = 1:length(ub)-1
    binVec = str2num(dec2bin(j)');
    binVec = [zeros(param.M-length(binVec),1);binVec];
    ub(:,j+1) = binVec .* masses;
end
param.ub = [ub,-ub(:,2:end)];
param.c = 1e-3*sqrt(sum(param.ub.^2,1));
param.gamma = 1e-2;

% generate divergence operator
param.div = sparse(param.edges(:),repmat(1:param.numE,1,2)',[ones(param.numE,1);-ones(param.numE,1)],param.N,param.numE);
param.gradDiv = param.div'*param.div;


%% solve transportation problem

u = ssn(param,@objfun_branchedTransport,zeros(param.M,param.numE));


%% visualize result

plot_flow(u,param);
