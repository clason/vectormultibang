function [J,G,H,mbInfo] = objfun_branchedTransport(param,x)
% OBJFUN_BRANCEHDTRANSPORT compute functional value, gradient, Hessian
% [J,G,H,AS] = OBJFUN_BRANCEHDTRANSPORT(PARAM,X) computes the value J of the functional to be
% minimized together with the gradient G and the Hessian H in the point X.
% MBINFO contains additional information on the multibang penalty and its
% regularization. The structure PARAM contains the problem parameters.
%
% X - M x N (N is number of streets, M number of materials)
% G - M x N
% H - operator from M x N arrays into M x N arrays
%
% April 12, 2021                    Christian Clason (c.clason@uni-graz.at)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

%% objective value
res = param.source - param.div * x';        % residual
J = res(:)'*res(:)/2 + sum(x.^2,1)*param.lengths*param.gamma/2 + param.auxOp.mb_penalty(x)*param.lengths;

%% gradient
[Hg,DHg,mbInfo.as,mbInfo.rnodes,~] = mb_general(res'*param.div./param.lengths',param.gamma,param.ub,param.c,param.auxOp);
G = x - Hg;

%% Hessian
H = @(dx) dx + permute(sum(DHg.*shiftdim(dx*param.div'*param.div./param.lengths',-1),2),[1 3 2]);
