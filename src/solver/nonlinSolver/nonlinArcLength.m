%=========================================================================%
% Fracture arc-length solver 
%=========================================================================%

function globdat = nonlinArcLength(props,globdat)

% Increment step and time
globdat.ts.step = globdat.ts.step + 1;

if ismember(globdat.ts.step,props.tStep.cutsteps) == 1
    globdat.ts.dt = globdat.ts.dt/props.tStep.cutsize;
end

globdat.ts.t    = globdat.ts.t + globdat.ts.dt;

% Display some information in command window
fprintf('\n')
disp('-------------------------------------------------------')
disp(['Step ', num2str(globdat.ts.step),', method: ',globdat.method])
disp(['nElems = ',num2str(globdat.mesh.numElements),...
      ', nDOFs = ',num2str(globdat.ndofs),...
      ', beta = ',num2str(globdat.beta)])

disp('-------------------------------------------------------')
disp(['Iteration    | ','Error'])

% Apply Constraints
globdat.state = applyConstraints(globdat.state,globdat.ts.t,globdat.consTable);

% Assemble system
globdat       = ModelAction(props,globdat,'AssembleCoupledSystem');

% Some more variables
globdat.Dstate = zeros(globdat.ndofs,1);
errVec         = zeros(props.nlSolver.maxiter,1);
iter           = 0;

% Enter iteration loop
while true

    % Increments and storage
    iter     = iter + 1;

    switch globdat.method

        case 'dispCtrl'

            delSol          = LinearSolver(globdat.CMat'*globdat.K*globdat.CMat, ...
                              globdat.CMat'*(-globdat.fint),props.linSolver);
            delSol          = globdat.CMat*delSol;

            globdat.state   = globdat.state  + delSol;
            globdat.Dstate  = globdat.Dstate + delSol;


        case 'arcLength'

            % Assemble augmented dissipation system
            globdat = ModelAction(props,globdat,'AssembleFractureDissipation');

            % Some temporary data (to be refined later)
            B1 = globdat.CMat'*globdat.K*globdat.hatSol;
            B2 = globdat.h'*globdat.CMat;
            K1 = globdat.CMat'*globdat.K*globdat.CMat; 
            
            f1 = globdat.CMat'*(-globdat.fint);
            f2 = -globdat.g + globdat.dTau;
        
            % Compute solution update
            dt   = (B2*(K1\f1)-f2)/(B2*(K1\B1));
            dSol = K1\(f1 - B1*dt);
            dSol = globdat.beta * globdat.CMat * dSol;
        
            % Update solution
            globdat.ts.dt = globdat.ts.dt + globdat.beta * dt;
            globdat.ts.t  = globdat.ts.t + globdat.beta * dt;
            globdat.state = globdat.state + dSol;
            globdat.state(globdat.loadDofs) = globdat.ts.t;

            globdat.Dstate  = globdat.state - globdat.state0;

        otherwise
            error('Bug in the code!')
    end

    % Assemble system
    globdat       = ModelAction(props,globdat,'AssembleCoupledSystem');

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
                                  
            % Do not accept the current load increment
            globdat.ts.step = globdat.ts.step - 1;

            % Set stepsize to zero only for arc-length
            if strcmp(globdat.method,'arcLength') == 1
                globdat.dt = 0;
            end

            % Refine and solution transfer
            globdat = refineNTransferSolution(props,globdat);

            % Recompute constraint matrix
            globdat = getConstraints(props,globdat);

            % Revert
            % globdat      = ModelAction(props,globdat,'revert');
            % globdat.ts.t = globdat.ts.t0;

            globdat.Dstate0 = globdat.Dstate;


            break;

        end

        %-------------------------------------------------------------%
        % If convergence is achieved and mesh is unchanged
        %-------------------------------------------------------------%

        % Compute dissipation
        globdat = ModelAction(props,globdat,'AssembleFractureDissipation');
        disp(['Dissipated Energy = ',num2str(globdat.g)])    
        
        % Re-size error vector and store in globdat
        errVec(iter+1:props.nlSolver.maxiter) = [];
        globdat.errVec{globdat.ts.step,1} = errVec;

        % Perform additional step depending of method
        switch globdat.method

            case 'dispCtrl'

                if globdat.g > props.nlSolver.arcLengthSwitchE
        
                    globdat.method = 'arcLength';
                    globdat.dTau   = props.nlSolver.arcLengthSwitchE;
                    globdat.ts.dt  = 0.0;

                    fprintf('\n')
                    disp('===============================')
                    disp('Switching to Arc-Length method!')
                    disp(['Setting dTau to ',num2str(globdat.dTau)])
                    disp('===============================')
                    fprintf('\n')
        
                end

            case 'arcLength'

                % Increase dTau and beta
                globdat.ts.dt = 0;
                if iter < props.nlSolver.arcLengthOptIter
                    globdat.dTau = 2.0 * globdat.dTau;
                    globdat.dTau = min(globdat.dTau,props.nlSolver.arcLengthMaxdTau);

                    disp(['Increasing dTau to ',num2str(globdat.dTau)])

                else
                    globdat.dTau = 0.5^(0.25*(iter-props.nlSolver.arcLengthOptIter));
                    globdat.dTau = min(globdat.dTau,props.nlSolver.arcLengthMaxdTau);

                    disp(['Setting dTau to ',num2str(globdat.dTau)])

                end
                globdat.beta = globdat.beta * 1.05;
                globdat.beta = min(globdat.beta,1.0);

            otherwise
                error('Bug in code!')
        end

        % Commit solution step
        globdat          = ModelAction(props,globdat,'commit');
        globdat.ts.t0    = globdat.ts.t;
        globdat.ts.fail  = 0;
        globdat.Dstate0  = globdat.Dstate;

        break;

        %-------------------------------------------------------------%
        % Terminate NR iterations when max iter is reached or when    
        % the error blows up                                          
        %-------------------------------------------------------------%
        elseif iter == props.nlSolver.maxiter || err > 5.0
            
            % Message to command window
            disp(['*** WARNING: Failed to converge in ',num2str(iter),...
                  ' iterations! ***'])

            % Do not accept the current step
            globdat.ts.step = globdat.ts.step - 1;

            % Perform additional step depending of method
            switch globdat.method
    
                case 'dispCtrl'

                    % Revert to old step and try with a smaller stepsize   
                    globdat.ts.t  = globdat.ts.t - globdat.ts.dt;
                    globdat.ts.dt   = globdat.ts.dt/props.tStep.cutsize;

                    disp(['Reducing stepsize to ',num2str(globdat.ts.dt)])
            
    
                case 'arcLength'
    
                    % Decrease dTau and beta
                    globdat.ts.dt = 0;

                    globdat.ts.fail = globdat.ts.fail + 1;

                    if globdat.ts.fail > 2

                        globdat.dTau = globdat.dTau/2.0; % 2.0 %dTau/1.25;
                        globdat.beta = globdat.beta/1.25;

                        disp(['Reducing dTau to ',num2str(globdat.dTau)])

                    end

                    globdat.beta = globdat.beta/1.25;
                    globdat.beta = max(globdat.beta,0.005);
    
                otherwise
                    error('Bug in code!')
            end

            % Undo all commit changes
            globdat          = ModelAction(props,globdat,'revert');
            globdat.ts.t     = globdat.ts.t0;
            globdat.Dstate   = globdat.Dstate0;

            break;

    % Exit convergence loop
    end

% Exit iteration loop    
end


end
%=========================================================================%
% End of function
%=========================================================================%

