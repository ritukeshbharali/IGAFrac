%=========================================================================%
% Setup problem (creates a globdat struct to store runtime variables)
%=========================================================================%

function globdat = setupProblem(props)

%-------------------------------------------------------------------------%
% Initialize globdat as struct
%-------------------------------------------------------------------------%
globdat = struct;


%-------------------------------------------------------------------------%
% Create mesh and compute Basis functions
%-------------------------------------------------------------------------%

[PHTelem,ctrlPts,mesh] = modelMesh(props.geom);
Basis                  = cartdev(PHTelem,ctrlPts,props.geom);

globdat.PHTelem        = PHTelem;
globdat.ctrlPts        = ctrlPts;
globdat.mesh           = mesh;
globdat.Basis          = Basis;


%-------------------------------------------------------------------------%
% Get FE Model
%-------------------------------------------------------------------------%
switch props.femodel.type

    case 'Elasticity'

        globdat.dim = props.geom.dim;
        assert(globdat.dim == length(props.femodel.dofs),'Elasticity FE Model: Dimension of problem does not match with dofs!')
        
        % Compute system dofs
        globdat.ndofs = globdat.dim*mesh.sizeBasis;


    case 'PhaseField'

        globdat.dim = props.geom.dim;
        assert(globdat.dim == length(props.femodel.dofs) - 1 ,'Phase-field FE Model: Dimension of problem does not match with dofs!')

        % Compute system dofs
        globdat.ndofs = (globdat.dim+1) * mesh.sizeBasis;

    otherwise
        error('FEModel: Not yet implemented!')
end


% Initialize state vectors (current, old step, old old step)
globdat.state        = zeros(globdat.ndofs,1);
globdat.state0       = zeros(globdat.ndofs,1);
globdat.state00      = zeros(globdat.ndofs,1);

% Initialize state increment vectors (current, old step)
globdat.Dstate       = zeros(globdat.ndofs,1);
globdat.Dstate0      = zeros(globdat.ndofs,1);

% Initialize int force vectors (current, old step)
globdat.fint         = zeros(globdat.ndofs,1);
globdat.fint0        = zeros(globdat.ndofs,1);


%-------------------------------------------------------------------------%
% Get and process model constraints
%-------------------------------------------------------------------------%

% Check if model has constraints
% ncons = length(props.cons);
% assert(ncons>0,...
%    'Without constraints, the system of equations would be singular!')

globdat = getConstraints(props,globdat);


%-------------------------------------------------------------------------%
% Initialize more variables
%-------------------------------------------------------------------------%

% Initialise vector for load-displacement plot
globdat.lodi            = zeros(props.tStep.nsteps+1,2);

% Parameters for timestepping and solution process
globdat.ts.step         = 0;
globdat.ts.t            = 0;
globdat.ts.dt           = props.tStep.initstepsize;

if strcmp(props.nlSolver.type,'qNonlin')
    globdat.ts.dt0      = props.tStep.initstepsize;
    globdat.ts.dt00     = props.tStep.initstepsize;
end

if strcmp(props.nlSolver.type,'nonlinArcLength')
    globdat.method      = 'dispCtrl';
    globdat.ts.t0       = 0.0;
    globdat.ts.fail     = 0;
end

% Line-search parameter
globdat.beta            = 1.0;

% Activate globdat for initiate timestepping
globdat.active          = true;
globdat.redo            = false;

% Create a cell to store error vectors for each step
globdat.errVec = cell(props.tStep.nsteps,1);

end

%=========================================================================%
% End of function
%=========================================================================%
