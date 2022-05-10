function globdat = getConstraints(props,globdat)

% Get the constraints table
try
    globdat.consTable = DirichletConstraints(globdat.PHTelem,props.geom,...
                                  globdat.dim*globdat.mesh.sizeBasis);
catch
    globdat.consTable = DirichletConstraints(globdat.PHTelem,props.geom,...
                                  globdat.dim*globdat.mesh.sizeBasis,...
                                  globdat.ctrlPts);
end

% Get load dofs
tmp                = find(globdat.consTable(:,2)==1);
globdat.loadDofs   = globdat.consTable(tmp,1);

% Build the constraint matrix

switch props.nlSolver.type

    case 'staggered'

        % Constraints for system 1
        tmp           = find(globdat.consTable(:,1)<=globdat.dim*globdat.mesh.sizeBasis);
        consT1        = globdat.consTable(tmp,:);

        % Constraint matrix for system 1
        globdat.CMat1 = speye(globdat.dim*globdat.mesh.sizeBasis);
        globdat.CMat1(:,consT1(:,1)) = [];

        % Constraints for system 1
        tmp           = find(globdat.consTable(:,1) > globdat.dim*globdat.mesh.sizeBasis);
        consT2        = globdat.consTable(tmp,:) - globdat.dim*globdat.mesh.sizeBasis;

        % Constraint matrix for system 2
        globdat.CMat2 = speye(globdat.mesh.sizeBasis);
        % globdat.CMat2(:,consT2(:,1)) = [];

    case 'acclStaggered'

        % Constraints for system 1
        tmp           = find(globdat.consTable(:,1)<=globdat.dim*globdat.mesh.sizeBasis);
        consT1        = globdat.consTable(tmp,:);

        % Constraint matrix for system 1
        globdat.CMat1 = speye(globdat.dim*globdat.mesh.sizeBasis);
        globdat.CMat1(:,consT1(:,1)) = [];

        % Constraints for system 1
        tmp           = find(globdat.consTable(:,1) > globdat.dim*globdat.mesh.sizeBasis);
        consT2        = globdat.consTable(tmp,:) - globdat.dim*globdat.mesh.sizeBasis;

        % Constraint matrix for system 2
        globdat.CMat2 = speye(globdat.mesh.sizeBasis);
        % globdat.CMat2(:,consT2(:,1)) = [];    

    case 'nonlinArcLength'

        globdat.CMat                           = speye(globdat.ndofs);
        globdat.CMat(:,globdat.consTable(:,1)) = [];

        % Hat solution vector
        globdat.hatSol                         = zeros(globdat.ndofs,1);
        globdat.hatSol(globdat.loadDofs)       = 1.0;

    otherwise

        globdat.CMat                           = speye(globdat.ndofs);
        globdat.CMat(:,globdat.consTable(:,1)) = [];

end

end