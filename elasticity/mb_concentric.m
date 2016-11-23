function [Hg,DHg,as,rnodes] = mb_concentric(p,alpha,gamma,ub)
% MB_CONCENTRIC compute regularized subdifferential, Newton derivative of multibang penalty (concentric corners)
% [HG,DHG,AS,RNODES] = MB_RADIAL(P,ALPHA,GAMMA,UB) computes the regularized 
% subdifferential HG, its Newton derivative DHG and the corresponding 
% active sets AS of the multibang penalty for concentric control states at 
% the dual variable P. RNODES is the number of nodes in active sets 
% corresponding to regularized singular arcs. ALPHA is the multibang 
% penalty parameter, GAMMA the Moreau-Yosida  regularization parameter and 
% UB the vector of admissible control states. 
%     
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

Q1 = p(:,1);  Q2 = p(:,2);  Qmag = sum(abs(p),2);  Qsup = max(abs(p),[],2); 
N2 = size(Q1,1);

%% compute active sets
eta = @(x) gamma*(x < 3*alpha + gamma) ...
    + (x-3*alpha).*(3*alpha + gamma <= x).*(x <= 3*alpha + 2*gamma)...
    + 2*gamma*(x > 3*alpha + 2*gamma);

iq = sign(Q1).*(abs(Q1) > eta(abs(Q2)));
jq = sign(Q2).*(abs(Q2) > eta(abs(Q1)));
kq = -1*(Qsup < 3*alpha +   gamma & Qmag < 3*alpha + 2*gamma) + ...
      1*(Qsup > 3*alpha + 2*gamma | Qmag > 3*alpha + 4*gamma);

as = zeros(N2,19);

% Q_ijk^gamma, i,j,k \neq 0
as(:,1) = (iq == 1)  & (jq == 1)  & (kq == -1);
as(:,2) = (iq == 1)  & (jq == -1) & (kq == -1);
as(:,3) = (iq == -1) & (jq == 1)  & (kq == -1);
as(:,4) = (iq == -1) & (jq == -1) & (kq == -1);

as(:,5) = (iq == 1)  & (jq == 1)  & (kq == 1);
as(:,6) = (iq == 1)  & (jq == -1) & (kq == 1);
as(:,7) = (iq == -1) & (jq == 1)  & (kq == 1);
as(:,8) = (iq == -1) & (jq == -1) & (kq == 1);

% Q_ijk^gamma, two zeros
as(:,9) = abs(iq) + abs(jq) + abs(kq) == 1;

% Q_ijk^gamma, one zero
as(:,10) = (iq ~= 0) & (jq == 0) & (kq ~= 0);
as(:,11) = (iq == 0) & (jq ~= 0) & (kq ~= 0);
as(:,12) = (iq ~= 0) & (jq ~= 0) & (kq == 0);

% number of nodes in regularized active sets
rnodes = nnz(as(:,9:end)); 

%% compute H_gamma(q)
Hg = as(:,1:8)*ub' + ...
     bsxfun(@times,(p - 3*alpha*[iq,jq])/gamma,as(:,9)) + ...
     bsxfun(@times,[(kq+3).*iq/2, Q2/gamma, ],as(:,10)) + ...
     bsxfun(@times,[Q1/gamma, (kq+3).*jq/2],as(:,11)) + ...
     bsxfun(@times,bsxfun(@times,(Qmag - 3*alpha),[iq,jq]/(2*gamma)),as(:,12));
        
%% compute DH_gamma(q)
DHg11 = sum(as(:,[9,11]),2)/gamma + as(:,12)/(2*gamma);
DHg12 = iq.*jq.*as(:,12)/(2*gamma);
DHg22 = sum(as(:,[9,10]),2)/gamma + as(:,12)/(2*gamma);

DHg = [spdiags(DHg11,0,N2,N2), spdiags(DHg12,0,N2,N2); ...
       spdiags(DHg12,0,N2,N2), spdiags(DHg22,0,N2,N2)];
end