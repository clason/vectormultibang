function [Hg,DHg,as,rnodes,auxOp] = mb_general(p,gamma,ub,c,auxOp)
% MB_GENERAL compute regularized subdifferential, Newton derivative of multibang penalty (general)
% [HG,DHG,AS,RNODES,AUXOP] = MB_GENERAL(P,GAMMA,UB,C,AUXOP) computes
% the regularized subdifferential HG, its Newton derivative DHG and the
% corresponding active sets AS of the multibang penalty for arbitrarily
% distributed control states at the dual variable P. RNODES is the number
% of nodes in active sets corresponding to regularized singular arcs.
% GAMMA is the Moreau-Yosida regularization parameter, UB the vector of
% admissible control states, C the vector of associated costs.
% AUXOP is a precomputed auxiliary operator; if empty, it will be recomputed.
% This implementation only covers the case where the graph of the multibang
% penalty solely consists of simplices.
%
% P - m x ...
% UB - m x M
% C - 1 x M
% Hg - m x ...
% DHg - m x m x ...
%
% April 19, 2021                    Christian Clason (c.clason@uni-graz.at)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

% set up matrices occurring inside Hg and DHg and store them in auxOp
m = size(ub,1);
if isempty(auxOp)
    auxOp = mb_general_preparation(gamma,ub,c,ndims(p)-1,[]);
end

% compute all possibilities for the prox operator (first m components)
% and the associated convex combination coefficients (last m+1 components)
permute1 = [auxOp.numDims+2 1 auxOp.numDims+3 2:auxOp.numDims+1];
permute2 = [1 4:3+auxOp.numDims 3 2];
wLambda = permute( sum( auxOp.A .* permute( p, permute1 ), 2 ) + auxOp.b, permute2 );

% compute active sets
tol = eps;
sz = size(p);
as = false([length(auxOp.regions) sz(2:auxOp.numDims+1)]);
for j = 1:length(auxOp.regions)
    as(j,auxOp.allDims{:}) = all( (wLambda(m+1:end,auxOp.allDims{:},j)>=0) & (wLambda(m+1:end,auxOp.allDims{:},j)<=1), 1 ) ...
        & (auxOp.gStar(wLambda(1:m,auxOp.allDims{:},j))-auxOp.gStarI(wLambda(1:m,auxOp.allDims{:},j),auxOp.regions{j}(1))<tol);
    as(j,auxOp.allDims{:}) = as(j,auxOp.allDims{:}) & ~any(as(1:j-1,auxOp.allDims{:}),1); % uses that auxOp.regions are sorted by number of elements
end
rnodes = nnz(as(1:auxOp.numReg,auxOp.allDims{:}));

% compute Hg
Hg = p;
for j = 1:length(auxOp.regions)
    aux = wLambda(1:m,auxOp.allDims{:},j);
    Hg(:,as(j,auxOp.allDims{:})) = Hg(:,as(j,auxOp.allDims{:})) - aux(:,as(j,auxOp.allDims{:}));
end
Hg = Hg/gamma;

% compute DHg
DHg = repmat(eye(m),[1 1 sz(2:auxOp.numDims+1)]);
for j = 1:length(auxOp.regions)
    DHg(:,:,as(j,auxOp.allDims{:})) = DHg(:,:,as(j,auxOp.allDims{:})) - auxOp.A(1:m,:,j);
end
DHg = DHg/gamma;

end
