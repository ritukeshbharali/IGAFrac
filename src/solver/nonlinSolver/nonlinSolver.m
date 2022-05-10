% Conventional Newton-Raphson solver 

function globdat = nonlinSolver(props,globdat)

% Increment step and time
globdat.ts.step = globdat.ts.step + 1;

if ismember(globdat.ts.step,props.tStep.cutsteps) == 1
    globdat.ts.dt = globdat.ts.dt/props.tStep.cutsize;
end

globdat.ts.t    = globdat.ts.t + globdat.ts.dt;

fprintf('\n')
disp('-------------------------------------------------------')
disp(['Step ', num2str(globdat.ts.step),', Method: Newton-Raphson'])
disp(['Elements = ',num2str(globdat.mesh.numElements),...
      ', DOFs = ',num2str(globdat.ndofs)])
disp('-------------------------------------------------------')
disp(['Iteration    | ','Error'])

errVec = zeros(props.nlSolver.maxiter,1);


% Apply Dirichlet Boundary Condition
globdat.state = ApplyDBC(globdat.state,globdat.ts.t,globdat.consTable);

% Set iteration count to zero and allocate vector to store error
iter   = 0;

% Enter staggered solution process
while true
    
    % Iteration counter increment, and store old internal force
        iter     = iter + 1;
        
        globdat         = ModelAction(props,globdat,'AssembleCoupledSystem');
        delSol          = LinearSolver(globdat.CMat'*globdat.K*globdat.CMat, ...
                          globdat.CMat'*(-globdat.fint),props.linSolver);
        delSol          = globdat.CMat*delSol;
        globdat.state   = globdat.state  + delSol;
                
        % Compute error and store in errCol vector
        if iter == 1
            err0 = norm(globdat.CMat'*globdat.fint);
            err = 1;
        else
            err = norm(globdat.CMat'*globdat.fint)/max(1,err0);
        end
        errVec(iter) = err;
        
        % Display current iteration and error info
        disp([num2str(iter),'            |',num2str(err)])
        
        % Check for convergence
        if err < props.nlSolver.tol
        
            % Display message to command window            
            fprintf('\n')
            disp(['Converged in ',num2str(iter),' iteration(s).'])
            loads = abs(sum(globdat.fint(globdat.loadDofs)));
            disp(['Load = ',num2str(loads),' Displacement = ',num2str(globdat.ts.t)])
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

        % Compute dissipation
        globdat = ModelAction(props,globdat,'AssembleFractureDissipation');
        disp(['Dissipated Energy = ',num2str(globdat.g)])    
            
        globdat = FEModelAction(props,globdat,'commit');
        
        % Re-size error vector and store in globdat
        errVec(iter+1:props.nlSolver.maxiter) = [];
        globdat.errVec{globdat.ts.step,1} = errVec;
        
        break;
            
        %-------------------------------------------------------------%
        % Terminate NR iterations when max iter is reached or when    %
        % the error blows up                                          %
        %-------------------------------------------------------------%
        elseif iter == props.nlSolver.maxiter || err > 1000.0
            
            % Message to command window
            disp('**** WARNING: Failed to converge! ****')
            
            % Revert variables to oldStep
            globdat.state  = globdat.state0;
            
            % Reduce the stepsize for prescribed displacement method
            globdat.ts.step = globdat.ts.step - 1;
            globdat.ts.t    = globdat.ts.t - globdat.ts.t;
            globdat.ts.dt   = globdat.ts.dt/props.tStep.cutsize;
                    
            break;
            
        end

end


end % eof

