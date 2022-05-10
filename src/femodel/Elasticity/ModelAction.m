%-------------------------------------------------------------------------%
% ModelAction is a wrapper function that carries out the action requested
% by the solver. Typically, these are assembling the matrices or internal
% force vectors, or swapping of history variables once converged is
% achieved.
%-------------------------------------------------------------------------%

function globdat = ModelAction(props,globdat,action)

switch action
        
    case 'AssembleCoupledSystem'
        
        [globdat.K, globdat.fint, globdat.markRef] = ...
            AssembleCoupledSystem(props,globdat);
        
    case 'AssembleLumpedMass'
        
        globdat.M = ...
            AssembleLumpedMass(props,globdat);
        
    case 'AssembleConsistentMass'
        
        globdat.M = ...
            AssembleConsistentMass(props,globdat);     
        
    case 'commit'    
        
        % Nothing to do
        

end

end



%-------------------------------------------------------------------------%
% Assemble coupled system                                                 %
%-------------------------------------------------------------------------%

function [K, fint, markRef] = AssembleCoupledSystem(props,globdat)

Sol      = globdat.state;
PHTelem  = globdat.PHTelem;
MeshInfo = globdat.MeshInfo;
Basis    = globdat.Basis;

geometry = props.geometry;
material = props.material;


% Some useful variables
dim        = geometry.dim;
nstress    = geometry.nstress;
numPatches = geometry.numPatches;
ngaussX    = geometry.ngaussX;
ngaussY    = geometry.ngaussY;
NN         = dim*MeshInfo.sizeBasis;

% Initialise cell for mesh refinement marker
markRef    = cell(numPatches,1);

% Initialise vector for global residual
fint      = zeros(NN,1);

% Initialise cells for local stiffness matrices
nelem      = MeshInfo.numElements;
kUU        = cell(nelem,1);

% Set index counter to zero
elementCounter = 0;
indexCounter   = 0;

% Loop over patches
for indexPatch = 1:numPatches
    
    % Initialise vectors for elements in a patch
    markRef(indexPatch) = {zeros(1,length(PHTelem{indexPatch}))};
    
    % Loop over elements
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            
            % Extract element displacement and phase-field
            nument       = size(PHTelem{indexPatch}(i).C,1);
            elemNodes    = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            elemUDofs    = reshape([2*elemNodes-1;2*elemNodes],1,dim*nument);
            elemDisp     = Sol(elemUDofs);
            
            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter + dim^2*nument^2;
            
            % Set local stiffness matrix, residual, displacement to zero
            localkUU      = zeros(dim*nument,dim*nument);
            localfintU    = zeros(dim*nument,1);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    % Get shape functions derivatives
                    [Bu,~,D]=strainGrad(Basis.dgdx{indexPatch}{i},nument,nstress,dim,kgauss,material.C);
                    
                    % Compute strain 
                    strain  = Bu*elemDisp;
                                        
                    % Compute stress                    
                    stress = material.C*strain;
                                            
                    % Compute local stiffness matrix
                    localkUU     = localkUU + (Bu'*D).* Basis.volume{indexPatch}{i}(kgauss);
                    
                    % Compute local fintU
                    localfintU    = localfintU + Bu' * stress .* Basis.volume{indexPatch}{i}(kgauss);
                    
                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            fint(elemUDofs)    = fint(elemUDofs) + localfintU;
            
            % Store element matrices in pre-allocated cells
            kUU{elementCounter}     = localkUU;
        end
    end
end
            
clear dsctrx

% COMMENT: Can we include the indexing in the previous loop?
II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S  = zeros(1,indexCounter);

% Set element and index counter to zero
elementCounter = 0;
indexCounter   = 0;

% TO-DO: Can we move this inside the previous loop? We will have to 
% pre-allocate II,JJ,S.
for indexPatch = 1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter+1;
            nument    = size(PHTelem{indexPatch}(i).C,1);
            sctrx     = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            dsctrx = reshape([2*sctrx-1;2*sctrx],1,dim*nument); 
            elmtStiff = [kUU{elementCounter}];
            II(indexCounter+1:indexCounter+dim^2*nument^2) = repmat(dsctrx,1,dim*nument);
            JJ(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(repmat(dsctrx,dim*nument,1),1,dim^2*nument^2);
            S(indexCounter+1:indexCounter+dim^2*nument^2)  = reshape(elmtStiff,1,dim^2*nument^2);
            indexCounter = indexCounter +dim^2*nument^2;
        end
    end
end

% Assemble (sparse) stiffness matrix
K = sparse(II,JJ,S,NN,NN);

end


%-------------------------------------------------------------------------%
% Assemble Lumped Mass                                                    %
%-------------------------------------------------------------------------%

function M = AssembleLumpedMass(props,globdat)

M       = [];

end




%-------------------------------------------------------------------------%
% Assemble Consistent Mass                                                %
%-------------------------------------------------------------------------%

function M = AssembleConsistentMass(props,globdat)

M       = [];

end



