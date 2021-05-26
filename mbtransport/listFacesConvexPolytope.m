function faces = listFacesConvexPolytope( A, b )
% lists all faces of the convex polytope defined by all x for which Ax<=b
% implements the recursive algorithm from
% Fukuda, Liebling, Margot: Analysis of backtrack algorithms for listing all vertices and all faces of a convex polyhedron
% see also https://www.cs.mcgill.ca/~fukuda/soft/polyfaq/node16.html
% A - m x n
% b - m x 1
% faces - cell array, containing lists of indices whose inequalities define a face, sorting decreasing by number of indices

% generate all faces
faces = faceEnum( A, b, [], [] );
% sort faces
lengths = zeros(length(faces),1);
for j = 1:length(faces)
    lengths(j) = length(faces{j});
end
[~,idx] = sort(lengths,1,'descend');
faces = faces(idx);
end

function exists = RFP( A, b, R, S )
[m,n] = size(A);
auxMat = [zeros(m,1) A];
auxMat(S,1) = 1;
Aeq = auxMat(R,:);
beq = b(R);
auxMat(R,:) = [];
b(R) = [];
options = optimoptions('linprog','Display','off');
[~,fval,exitflag,~] = linprog([-1;zeros(n,1)],auxMat,b,Aeq,beq,[],[1;Inf*ones(n,1)],options);
if isempty(fval) || ( exitflag <= 0 )
    exists = false;
else
    exists = ( fval < 0 );
end
end

function face = minimalFace( A, b, R, S )
face = R;
J = (1:size(A,1))';
J(R) = 0;
J(S) = 0;
J(J==0) = [];
for j = J
    if ~RFP( A, b, R, j )
        face = [face;j];
    end
end
end

function faces = faceEnum( A, b, R, S )
if RFP( A, b, R, S )
    F = minimalFace( A, b, R, S );
    J = (1:size(A,1))';
    J(F) = 0;
    J(S) = 0;
    J(J==0) = [];
    faces = {F};
    for k = 1:length(J)
        faces = [faces,faceEnum( A, b, [F;J(k)], [S;J(1:k-1)] )];
    end
else
    faces = {};
end
end
