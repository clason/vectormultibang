function x = ssn(param,objfun,x0)
% SSN semi-smooth Newton method with path-following and line search
% X = SSN(PARAM,OBJFUN,X0) computes the optimal material flux using a
% semi-smooth method. The functional to be minimized is specified using the
% function handle OBJFUN, which evaluates functional, gradient and Hessian.
% The structure PARAM contains the problem parameters. X0 is the initial
% guess for the control X.
%
% April 12, 2021                    Christian Clason (c.clason@uni-graz.at)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

maxit  = 20;                      % maximum number of SSN iterations
reltol = 1e-9;                    % relative tolerance for gradient norm
cgtol  = 1e-11;                   % desired reduction of residual in CG
cgits  = param.numE;              % maximum number of CG iterations
x = x0;                           % initial guess

warning('off','MATLAB:gmres:tooSmallTolerance');
%% continuation in gamma
param.gamma = 2e1;                % initial Moreau--Yosida parameter
q = 0.5;                          % reduction parameter
param.auxOp.regions = [];
while param.gamma > 1e-7
    fprintf('\nCompute solution for gamma = %1.3e (1-q = %1.3e):\n',param.gamma,1-q);
    
    % update operators for regularized multibang derivatives
    param.auxOp = mb_general_preparation(param.gamma,param.ub,param.c,ndims(x0)-1,param.auxOp.regions);
    
    fprintf('Iter    objective | normgrad    dAS | stepsize   relres  CGit\n');
    
    it = 0;  GGold = 1e99;  abstol = min(1e-6,param.gamma);
    as_old = zeros(length(param.auxOp.regions),param.numE);  tau = 1;  dx = zeros(size(x));
    while it <= maxit
        % compute new gradient
        [J,G,H,mbInfo] = objfun(param,x);
        if it == 0
            G0 = sqrt(tau)*norm(G);
            flag = 0;  cgit = 0;
        end
        
        % line search on gradient norm (correctly scaled discrete norm)
        GG = sqrt(tau)*norm(G);
        if GG >= GGold       % if no decrease: backtrack (never on iteration 1)
            tau = tau/2;
            x = x - tau*dx;
            if tau >= 1e-5   % bypass rest of loop; backtrack further
                continue;
            end
        end
        
        % compute statistics and change in active sets
        as_change = nnz(mbInfo.as~=as_old)/2; % number of points that changed
        
        % output iteration details
        fprintf('%3d:  %1.5e | %1.3e  ', it, J, GG);
        if it > 0
            fprintf('%4d | %1.1e   %1.1e   %d\n', as_change, tau, relres, cgit(2));
        else
            fprintf('\n');
        end
        
        % terminate Newton?
        if (GG < reltol*sqrt(G0)) && (as_change == 0)  % convergence (relative norm)
            fprintf('\n#### converged with relative tol: |grad|<=%1.1e |grad0|\n',reltol);
            flag = 0;
            break;
        elseif (GG < abstol) && (as_change == 0)  % convergence (absolute norm)
            fprintf('\n#### converged with absolute tol: |grad|<=%1.1e\n',abstol);
            flag = 1;
            break;
        elseif it == maxit                     % failure, too many iterations
            fprintf('\n#### not converged: too many iterations\n');
            flag = 2;
            break;
        elseif tau < 1e-5                      % failure, too small stepsize
            fprintf('\n#### not converged: too small stepsize\n');
            flag = 3;
            break;
        end
        
        % otherwise update information, continue
        it = it+1;  GGold = GG;  as_old = mbInfo.as;  tau = 1;
        % compute Newton step, update
        hessianOp = @(dx) reshape(H(reshape(dx,param.M,[])),[],1);
        [dx, flag, relres, cgit] = gmres(hessianOp, -G(:), [], cgtol, cgits);
        dx = reshape(dx,param.M,[]);
        x = x + dx;
    end
    
    if flag >= 2
        param.gamma = param.gamma/q;
        q = q^0.25;
        param.gamma = param.gamma*q;
        if 1-q > 1e-6
            x = x0;
        end
    else
        if it <= maxit*1/4
            q = min(1-1e-3,max(q^1.25,0.5));
        elseif it >= maxit*3/4
            q = min(1-1e-4,q^0.75);
        else
            q = min(1-1e-3,q);
        end
        param.gamma = param.gamma*q;
        x0 = x;
    end
    
    % output number of nodes in regularized active sets
    fprintf('Solution has %i node(s) in regularized active sets\n', mbInfo.rnodes);
    
end

