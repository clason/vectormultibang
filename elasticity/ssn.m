function [u,y] = ssn(z,A,M,N,mb_penalty,uplot,yplot)
% SSN semi-smooth Newton method with continuation and line search
% [U,Y] = SSN(Z,A,M,N,mb_penalty,pplot) computes the optimal control U and the 
% corresponding state Y using a semi-smooth method. 
%
% November 21, 2016          Christian Clason (christian.clason@uni-due.de)
%              Carla Tameling (carla.tameling@mathematik.uni-goettingen.de)
%                           Benedikt Wirth (benedikt.wirth@uni-muenster.de)

%% setup 
maxit = 50;               % max number of Newton steps
% precompute some terms
MTz = M*z(:);    AT = A';    N2 = N*N;
% initialize iterates
y  = zeros(N2,2);                      % state variable
p  = zeros(N2,2);                      % dual variable
[~,~,as,~] = mb_penalty(p,1); as_old = 0*as; % nasty hack to get right shape
%% continuation in gamma
gamma  = 1e2;
while gamma > 1e-10
    it = 1;    nold = 1e99;    tau = 1;
    
    fprintf('\nCompute solution for gamma = %1.3e:\n',gamma);
    fprintf('Iter  |  normgrad      dAS     stepsize\n');
    while true
        % update active sets 
        [Hg,DHg,as,rnodes] = mb_penalty(p,gamma);
 
        % right hand side
        C   = [M AT; A -M*DHg];
        rhs = [MTz-M*y(:)-AT*p(:); -A*y(:) + M*Hg(:)];
        nr  = norm(rhs(:));
        
        % line search
        if nr >= nold             % if no decrease: backtrack
            tau = tau/2;
            y(:) = y(:) - tau*dx(1:2*N2);
            p(:) = p(:) - tau*dx(1+2*N2:end);
            if tau < 1e-6         % accept non-monotone step
            else                  % bypass rest of while loop;
                continue;
            end
        end
        
        % terminate Newton?
        update = nnz((as-as_old));
        fprintf('%2d    |  %1.3e   %5d     %1.2e\n',it,nr,update,tau);
        if update == 0  && nr < 1e-6  % success, solution found
            break;
        elseif it == maxit            % failure, too many iterations
            break;
        end        
        
        % semismooth Newton step
        dx = C\rhs;
        y(:) = y(:)+dx(1:2*N2);
        p(:) = p(:)+dx(2*N2+1:end);
        
        % otherwise update information, continue
        it = it+1;   nold = nr;   tau = 1;   as_old = as;
    end %newton
    
    % check convergence
    if it < maxit                      % converged: accept iterate
        u = Hg;
                
        fprintf('Solution has %i node(s) in regularized active sets\n',rnodes);
        
        if rnodes == 0 || it == 1      % solution optimal: terminate
            break;
        else                           % reduce gamma, continue
            gamma = gamma/2;
        end
    else                               % not converged: reject, terminate
        fprintf('Iterate rejected, returning u_gamma for gamma = %1.3e\n',...
            gamma*2);
        break;
    end

    % show current control and corresponding deformation
    uplot(u);  
    yplot(y);  
    drawnow update
end
