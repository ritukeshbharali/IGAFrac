%=========================================================================%
% Anderson accelerated alternate minization (staggered) solver
%=========================================================================%

function globdat = acclStaggeredSolver(props,globdat)

%-------------------------------------------------------------------------%
% Increment step and (pseudo) time
%-------------------------------------------------------------------------%

globdat.ts.step = globdat.ts.step + 1;

if ismember(globdat.ts.step,props.tStep.cutsteps) == 1
    globdat.ts.dt = globdat.ts.dt/props.tStep.cutsize;
end

globdat.ts.t    = globdat.ts.t + globdat.ts.dt;

%-------------------------------------------------------------------------%
% Print information to the command window
%-------------------------------------------------------------------------%

fprintf('\n')
disp('-------------------------------------------------------')
disp(['Step ', num2str(globdat.ts.step)])
disp(['Elements = ',num2str(globdat.mesh.numElements),...
      ', DOFs = ',num2str(globdat.ndofs)])
disp('-------------------------------------------------------')
disp(['Iteration    | ','Error'])

%-------------------------------------------------------------------------%
% Get some useful variables
%-------------------------------------------------------------------------%

nDof1    = globdat.dim*globdat.mesh.sizeBasis;
nDof2    = globdat.mesh.sizeBasis;
errVec   = zeros(props.nlSolver.maxiter,1);

dSol2    = zeros(nDof2,1);
olddSol2 = zeros(nDof2,1);

%-------------------------------------------------------------------------%
% Solve a step
%-------------------------------------------------------------------------%

% Apply Dirichlet Boundary Condition
globdat.state = ApplyDBC(globdat.state,globdat.ts.t,globdat.consTable);

% Set iteration count to zero and allocate vector to store error
iter   = 0;

% Anderson acceleration arrays
m    = 5;
mk   = 0;
beta = 1.0;
FkPf = [];
EkPf = [];

% Begin iterations
while true
    
    % Iteration counter increment, and store old internal force
    iter     = iter + 1;
    oldstate = globdat.state;
    olddSol2 = dSol2;
        
    % Solve system 1
    globdat                = ModelAction(props,globdat,'AssembleSystem1');
    if iter == 1
        dSol1                  = LinearSolver(globdat.CMat1'*globdat.K1*globdat.CMat1, ...
                             globdat.CMat1'*(-globdat.fint1),props.linSolver);
    else
        dSol1                  = LinearSolver(globdat.CMat1'*globdat.K10*globdat.CMat1, ...
                             globdat.CMat1'*(-globdat.fint1),props.linSolver);
    end
    dSol1                  = globdat.CMat1*dSol1;
    globdat.state(1:nDof1) = globdat.state(1:nDof1)  + dSol1;
            
    % Solve system 2
    globdat                = ModelAction(props,globdat,'AssembleSystem2');
    if length(globdat.CMat2(1,:)) < globdat.mesh.sizeBasis
        if iter == 1
            dSol2             = LinearSolver(globdat.CMat2'*globdat.K2*globdat.CMat2, ...
                                 globdat.CMat2'*(-globdat.fint2),props.linSolver);
        else
            dSol2             = LinearSolver(globdat.CMat2'*globdat.K20*globdat.CMat2, ...
                                 globdat.CMat2'*(-globdat.fint2),props.linSolver);
        end
        dSol2                  = globdat.CMat2*dSol2;
    else
        if globdat.ts.step == 1 && iter == 1
            dSol2              = LinearSolver(globdat.K2,-globdat.fint2,props.linSolver);
        else
            dSol2              = LinearSolver(globdat.K20,-globdat.fint2,props.linSolver);
        end
    end

    % Update based on Anderson acceleration
    if iter < 3 
        globdat.state(nDof1+1:nDof1+nDof2) = globdat.state(nDof1+1:nDof1+nDof2) + dSol2;
    else
        mk = mk+1;
            
            if mk <= m
            
                EkPf(:,mk) = globdat.state(nDof1+1:nDof1+nDof2) - oldstate(nDof1+1:nDof1+nDof2);
                FkPf(:,mk) = dSol2 - olddSol2;
            
            else
                
                EkPf = [EkPf,globdat.state(nDof1+1:nDof1+nDof2) - oldstate(nDof1+1:nDof1+nDof2)];
                FkPf = [FkPf,dSol2 - olddSol2];
                
                EkPf(:,1) = [];
                FkPf(:,1) = [];
                
                mk = m;
                
            end
            
            gamk = lsqlin(FkPf(:,mk),dSol2,[],[]);
            
            globdat.state(nDof1+1:nDof1+nDof2)  = globdat.state(nDof1+1:nDof1+nDof2) ...
                            + beta*dSol2 ...
                            -(EkPf(:,mk) + beta*FkPf(:,mk))*gamk ;
    end

    if iter == 1
        globdat.K10 = globdat.K1;
        globdat.K20 = globdat.K2;
    end
    
    % Compute error and store in errCol vector
    if iter == 1
        err0 = norm(globdat.state);
        err = 1;
    else
        err = norm(globdat.state-oldstate)/max(1,err0);
    end
    errVec(iter) = err;
    
    % Display current iteration and error info
    disp([num2str(iter),'            |',num2str(err)])
    
    %---------------------------------------------------------------------%
    % Check for convergence
    %---------------------------------------------------------------------%
    if err < props.nlSolver.tol

        % Display message to command window            
        fprintf('\n')
        disp(['Converged in ',num2str(iter),' iteration(s).'])

        loads = abs(sum(globdat.fint1(globdat.loadDofs)));
        disp(['Load = ',num2str(loads),' Displacement = ', ...
              num2str(globdat.ts.t)])
        
        % Store load-displacement data
        globdat.lodi(globdat.ts.step+1,:) = [loads,globdat.ts.t];

        %-------------------------------------------------------------%
        % Mesh refinement and solution transfer                       %
        %-------------------------------------------------------------%
            
        % Check if mesh needs to be refined
        if sum([globdat.markRef{:}])
                           
            % Display message to command window
            fprintf('\n')
            disp('**** MESH CHANGED, re-do step! ****')

            %-------------------------------------------------------------%
            % Refine and transfer nodal data
            %-------------------------------------------------------------%
            
            % Refine and solution transfer
            globdat = refineNTransferSolution(props,globdat);

            % Recompute constraint matrix
            globdat = getConstraints(props,globdat);
                                  
            % Do not accept the current load increment
            globdat.ts.step = globdat.ts.step - 1;
            globdat.ts.t    = globdat.ts.t - globdat.ts.dt;
            
            break;

        end
            
        globdat = FEModelAction(props,globdat,'commit');
        
        % Re-size error vector and store in globdat
        errVec(iter+1:props.nlSolver.maxiter) = [];
        globdat.errVec{globdat.ts.step,1} = errVec;
        
        break;
            
        %-------------------------------------------------------------%
        % Terminate NR iterations when max iter is reached or when    
        % the error blows up                                          
        %-------------------------------------------------------------%
        elseif iter == props.nlSolver.maxiter || err > 5.0
            
            % Message to command window
            disp(['*** WARNING: Failed to converge in ',num2str(iter),...
                  ' iterations! ***'])
            disp(['*** Reducing stepsize by factor ', ...
                  num2str(props.tStep.cutsize)])
            
            % Revert variables to oldStep
            globdat.state  = globdat.state0;
            
            % Reduce the stepsize for prescribed displacement method
            globdat.ts.step = globdat.ts.step - 1;
            globdat.ts.t    = globdat.ts.t - globdat.ts.dt;
            globdat.ts.dt   = globdat.ts.dt/props.tStep.cutsize;
                    
            break;
            
    end

end

end 

%=========================================================================%
% End of function
%=========================================================================%

