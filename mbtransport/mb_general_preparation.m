function auxOp = mb_general_preparation(gamma,ub,c,numDims,regions)
% MB_GENERAL_PREPARATION compute auxiliary operator encoding information on the structure of multibang penalty
% AUXOP = MB_GENERAL_PREPARATION(GAMMA,UB,C,REGIONS) computes an auxiliary
% operator AUXOP (a struct) for the multibang penalty (general).
% GAMMA is the Moreau-Yosida regularization parameter, UB the vector of
% admissible control states, C the vector of associated costs,
% NUMDIMS the number of spatial dimensions of the control variable to be optimized later,
% REGIONS a list of different domains of the multibang penalty or its conjugate;
% if empty, it is recomputed.
%
% UB - m x M
% C - 1 x M
%
% April 19, 2021                    Christian Clason (c.clason@uni-graz.at)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

% compute topological information of epigraph of multibang penalty
m = size(ub,1);
if isempty(regions)
    % list all faces of the epigraph of convex conjugate of multibang penalty; a face is described by the corresponding list of inequality indices
    % equivalently: list all faces of the epigraph of multibang penalty; a face is described by the corresponding list of vertices
    auxOp.regions = listFacesConvexPolytope([ub' -ones(length(c),1)],c');
    % remove full epigraph (the highest-dimensional face, corresponding to an empty index set)
    auxOp.regions = auxOp.regions(1:end-1);
else
    auxOp.regions = regions;
end
auxOp.numReg = length(auxOp.regions);   % number of the regularized faces of g*
while length(auxOp.regions{auxOp.numReg}) == 1
    auxOp.numReg = auxOp.numReg - 1;
end
auxOp.numFacets = 1;     % number of facets of the graph of g
while length(auxOp.regions{auxOp.numFacets+1}) > m
    auxOp.numFacets = auxOp.numFacets + 1;
end

% calculate matrices defining Hg and DHg
auxOp.A = zeros(2*m+1,m,length(auxOp.regions));
auxOp.b = zeros(2*m+1,1,length(auxOp.regions));
for j = 1:length(auxOp.regions)
    numIdx = length(auxOp.regions{j});
    mat = [eye(m) gamma*ub(:,auxOp.regions{j});
        (ub(:,auxOp.regions{j}(2:end))-ub(:,auxOp.regions{j}(1)))' zeros(numIdx-1,numIdx);
        zeros(1,m) ones(1,numIdx)];
    matInv = mat\eye(m+numIdx);
    auxOp.A(1:m+numIdx,:,j) = matInv(:,1:m);
    auxOp.b(1:m+numIdx,:,j) = matInv(:,m+1:end)*[c(auxOp.regions{j}(2:end))-c(auxOp.regions{j}(1)),1]';
end

% compute different functions (g,g*)
auxOp.numDims = numDims;
dimPermute = [1 3:2+auxOp.numDims 2];
ub = permute(ub,dimPermute);
c = permute(c,dimPermute);
auxOp.gStar = @(q) max( sum( ub .* q, 1 ) - c, [], auxOp.numDims+2 );
auxOp.gStarI = @(q,i) sum( ub(:,i) .* q, 1 ) - c(i);
auxOp.allDims = repmat({':'},1,auxOp.numDims);
auxOp.verticesGStar = auxOp.b(1:m,auxOp.allDims{:},1:auxOp.numFacets);
auxOp.vertexValuesGStar = shiftdim( zeros(1,auxOp.numFacets), -auxOp.numDims );
auxOp.getBarycentricCoordsOp = zeros( m+1, m+1, auxOp.numFacets );
for j = 1:auxOp.numFacets
    auxOp.vertexValuesGStar(j) = auxOp.gStarI(auxOp.verticesGStar(:,j),auxOp.regions{j}(1));
    mat = [ub(:,auxOp.regions{j});ones(1,m) 0];
    auxOp.getBarycentricCoordsOp(:,:,j) = mat\eye(m+1);
end
auxOp.getBarycentricCoordsOp = permute( auxOp.getBarycentricCoordsOp, [1 2 4:3+auxOp.numDims 3] );
inUnitInterval = @(f) (f>=0)&(f<=1);
auxOp.g = @(u) true2inf( shiftdim( all(any(~inUnitInterval(sum(auxOp.getBarycentricCoordsOp(:,1:m,auxOp.allDims{:},:).*shiftdim(u,-1),2)+auxOp.getBarycentricCoordsOp(:,m+1,auxOp.allDims{:},:)),1),auxOp.numDims + 3), 1 ) ) ...
    + max( sum( u .* auxOp.verticesGStar, 1 ) - auxOp.vertexValuesGStar, [], auxOp.numDims + 2 );
auxOp.mb_penalty = auxOp.g;
auxOp.verticesGStar = squeeze(auxOp.verticesGStar);
auxOp.vertexValuesGStar = shiftdim(squeeze(auxOp.vertexValuesGStar),-1);

end


function result = true2inf(f)
result = zeros(size(f));
result(f) = Inf;
end
