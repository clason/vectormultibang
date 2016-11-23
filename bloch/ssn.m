function x = ssn(d,objfun,applyHess)
% SSN semi-smooth Newton method with continuation and line search
% X = SSN(D,OBJFUN,APPLYHESS,X0) computes the optimal control X using a
% semi-smooth method. The functional to be minimized is specified using the
% function handles OBJFUN, which evaluates functional and gradient, and
% APPLYHESS, which computes the action of the Hessian on a given direction.
% The structure PARAM contains the problem parameters. X0 is the initial
% guess for the control X.
%
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

maxit  = 500;                     % maximum number of SSN iterations
reltol = 1e-7;                    % relative tolerance for gradient norm
abstol = 1e-7;                    % absolute tolerance for gradient norm
cgtol  = 1e-10;                   % desired reduction of residual in CG
cgits  = d.Nu;                    % maximum number of CG iterations
x = [d.u0;d.v0];                  % initial guess

%% continuation in gamma
d.gamma  = 1e2;                   % initial Moreau--Yosida parameter
while d.gamma > 1e-10
    
    fprintf('\nCompute solution for gamma = %1.3e:\n',d.gamma);
    fprintf('Iter    objective | normgrad    dAS | stepsize flag relres CGit\n');
    
    it = 0;  GGold = 1e99;
    as_old = zeros(d.Nu,4*d.M + 1);  tau = 1;  dx = zeros(size(x));
    while it <= maxit
        % compute new gradient
        [J,G,Xk] = objfun(d,x);
        if it == 0
            G0 = sqrt(tau)*norm(G);
            flag = 0;  cgit = 0;
        end
        
        % line search on gradient norm (correctly scaled discrete norm)
        GG = sqrt(tau)*norm(G);
        if GG >= GGold       % if no decrease: backtrack (never on iteration 1)
            tau = tau/2;
            x = x - tau*dx;
            if tau < 1e-6    % if step too small: accept nonmonotone step
                flag = 3;
            else             % else: bypass rest of loop; backtrack further
                continue;
            end
        end
        
        % compute statistics and change in active sets
        as_change = nnz(abs(Xk.As - as_old)>0.5); % number of points that changed
        
        % output iteration details
        fprintf('%3d:  %1.5e | %1.3e  ', it, J,  GG);
        if it > 0
            fprintf('%4d | %1.1e   %d  %1.1e    %d\n', as_change, tau, flag, relres, cgit(2));
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
        end
        
        % otherwise update information, continue
        it = it+1;  GGold = GG;  as_old = Xk.As;  tau = 1;
        % compute Newton step, update
        DG = @(dx) applyHess(d,Xk,dx);  % Hessian
        % compute Newton step, update
        [dx, flag, relres, cgit] = gmres(DG, -G, cgits, cgtol, cgits);
        x = x + dx;
    end
    
    % output number of nodes in regularized active sets
    fprintf('Solution has %i node(s) in regularized active sets\n', Xk.rnodes);
    
    % show current control, state
    blochplot(d,x);  drawnow update
    
    % reduce gamma
    d.gamma = d.gamma/2;
end

