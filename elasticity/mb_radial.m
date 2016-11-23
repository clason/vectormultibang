function [Hg,DHg,as,rnodes] = mb_radial(p,alpha,gamma,ub,ubins,omega0)
% MB_RADIAL compute regularized subdifferential, Newton derivative of multibang penalty (radial)
% [HG,DHG,AS,RNODES] = MB_RADIAL(P,ALPHA,GAMMA,UB,UBINS,OMEGA0) computes
% the regularized subdifferential HG, its Newton derivative DHG and the
% corresponding active sets AS of the multibang penalty for radially 
% distributed control states at the dual variable P. RNODES is the number 
% of nodes in active sets corresponding to regularized singular arcs.
% ALPHA is the multibang penalty parameter, GAMMA the Moreau-Yosida 
% regularization parameter, UB the vector of admissible control states. 
% UBINS are the boundaries of the corresponding (radial) Voronoi cells and 
% OMEGA0 is the fixed magnitude of the control states.
%     
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

Q1 = p(:,1);   Q2 = p(:,2);
M = size(ub,2); N2 = size(Q1,1); 
om02 = omega0^2;

cplus = @(i) mod(i,M)+1;                    % i+1 mod M
vprod = @(u) sum(bsxfun(@times,p,u),2);     % <q,ub_i> ub_i

%% compute active sets
Theta    = angle(Q1+1j*Q2);
[~,~,iq] = histcounts(Theta,ubins);   iq(iq==M+1)=1; % wrap around

Theta    = angle(Q1-gamma*ub(1,iq)' + 1j*(Q2-gamma*ub(2,iq)'));
[~,~,jq] = histcounts(Theta,ubins);   jq(jq==M+1)=1; % wrap around

uijq = ub(:,iq) + ub(:,jq);
sigq = (Q1-gamma/2*uijq(1,:)').*uijq(1,:)' + (Q2-gamma/2*uijq(2,:)').*uijq(2,:)';
rhoq = (Q1.*ub(1,iq)' + Q2.*ub(2,iq)');
qfac = rhoq/om02 - alpha/2;

Theta    = angle((Q1-qfac.*ub(1,iq)') + 1j*(Q2-qfac.*ub(2,iq)'));
[~,~,kq] = histcounts(Theta,ubins);   kq(kq==M+1)=1; % wrap around

as = zeros(N2,4*M+1);
for i = 1:M
    % Q_i gamma i > 0
    as(:,i) =   rhoq > (alpha/2+gamma)*om02 & iq==i & jq==i;
    % Q_i0 gamma
    as(:,i+M) = rhoq >= alpha/2*om02 & rhoq <= (alpha/2+gamma)*om02 & iq==i & kq==i;
    % Q_i,i+1 gamma
    as(:,i+2*M) = sigq >  alpha*om02 & ((iq==i & jq==cplus(i)) | (jq==i & iq==cplus(i)));
    % Q_0,i,i+1 gamma
    as(:,i+3*M) = sigq <= alpha*om02 & ((iq==i & kq==cplus(i)) | (kq==i & iq==cplus(i)));
end
% Q_0 gamma: all other nodes
as(:,4*M+1) = 1-sum(as(:,1:4*M),2);

% number of nodes in regularized active sets
rnodes = nnz(as(:,M+1:end-1)); 

%% compute H_gamma(q)
% Q_i^gamma
Hg = as(:,1:M)*ub';
for i = 1:M
    usum  = ub(:,i) + ub(:,cplus(i));   usumn = usum/sum(usum.^2); 
    udif  = ub(:,i) - ub(:,cplus(i));   udifn = udif/sum(udif.^2);
    % Q_0,i^gamma
    Hg = Hg + (vprod(ub(:,i)')/om02 - alpha/2).*as(:,M+i)*ub(:,i)'/gamma;
    % Q_i,i+1^gamma 
    Hg = Hg + bsxfun(@times,[usum(1)/2 + vprod(udif')*udifn(1)/gamma, ...
                             usum(2)/2 + vprod(udif')*udifn(2)/gamma], ...
                            as(:,i+2*M));
    % Q_0,i,i+1^gamma
    Hg = Hg + bsxfun(@times,[Q1 - alpha*om02*usumn(1), ...
                             Q2 - alpha*om02*usumn(2)]/gamma, ...
                            as(:,i+3*M));
end

%% compute DNH_gamma(q)
% Q_0,i gamma and Q_0,i,i+1 gamma
fac1 = 1/(gamma*om02)*as(:,(M+1):2*M);
diag11 = fac1*(ub(1,:).^2)' + 1/gamma*sum(as(:,(3*M+1):4*M),2);
diag12 = fac1*(ub(1,:).*ub(2,:))';
diag22 = fac1*(ub(2,:).^2)' + 1/gamma*sum(as(:,(3*M+1):4*M),2);

for i = 1:M % Q_i,i+1 gamma
    udif  = ub(:,i) - ub(:,cplus(i));
    fac2 = 1/(gamma * sum(udif.^2));
    diag11 = diag11 + fac2*as(:,i+2*M)*udif(1)^2;
    diag12 = diag12 + fac2*as(:,i+2*M)*udif(1)*udif(2);
    diag22 = diag22 + fac2*as(:,i+2*M)*udif(2)^2;
end

DHg11 = spdiags(diag11,0,N2,N2);
DHg12 = spdiags(diag12,0,N2,N2);
DHg22 = spdiags(diag22,0,N2,N2);

DHg = [DHg11 DHg12; DHg12 DHg22];
end
