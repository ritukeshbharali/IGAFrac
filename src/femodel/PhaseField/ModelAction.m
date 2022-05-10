%-------------------------------------------------------------------------%
% ModelAction is a wrapper function that carries out the action requested
% by the solver. Typically, these are assembling the matrices or internal
% force vectors, or swapping of history variables once converged is
% achieved.
%-------------------------------------------------------------------------%

function globdat = ModelAction(props,globdat,action)

switch action
    
    case 'AssembleSystem1'
        
        [globdat.K1, globdat.fint1, globdat.markRef] = ...
            AssembleDisplacementSystem(props,globdat);
        
    case 'AssembleSystem2'
        
        [globdat.K2, globdat.fint2] = ...
            AssemblePhaseFieldSystem(props,globdat);
        
    case 'AssembleCoupledSystem'
        
        [globdat.K, globdat.fint, globdat.markRef] = ...
            AssembleCoupledSystem(props,globdat);
        
    case 'AssembleFractureDissipation'
        
        [globdat.h, globdat.g] = ...
            AssembleFractureDissipation(props,globdat); 
        
    case 'AssembleLumpedMass'
        
        globdat.M = ...
            AssembleLumpedMass(props,globdat);
        
    case 'AssembleConsistentMass'
        
        globdat.M = ...
            AssembleConsistentMass(props,globdat); 
        
    case 'AssemblePhaseFieldMass'
        
        globdat.M = ...
            AssemblePhaseFieldMass(props,globdat);    
        
    case 'commit'    
        
        globdat.state00 = globdat.state0;
        globdat.state0  = globdat.state;

    case 'revert'    
        
        globdat.state   = globdat.state0;
        globdat.state0  = globdat.state00;    
        

end

end



%-------------------------------------------------------------------------%
% Assemble displacement system                                            %
%-------------------------------------------------------------------------%

function [KU, fintU, markRef] = AssembleDisplacementSystem(props,globdat)

Sol      = globdat.state;
PHTelem  = globdat.PHTelem;
MeshInfo = globdat.mesh;
Basis    = globdat.Basis;

geometry = props.geom;
material = props.mat;


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
fintU      = zeros(NN,1);

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
            elemPhiDofs  = elemNodes+NN;
            elemDisp     = Sol(elemUDofs);
            elemPhi      = Sol(elemPhiDofs);
            
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
                    
                    % Get phase-field and compute degradation function
                    phigp   = Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi;
                    
                    fac1   = (1-phigp)^material.p;
                    fac2   = fac1 + material.a1*phigp + material.a1*material.a2*phigp*phigp + ...
                               material.a1*material.a2*material.a3*phigp*phigp*phigp;
                    gphi   = fac1/fac2;                    
                    
                    % Set refinement based on phase-field
                    if (phigp > geometry.threshPhi) && (PHTelem{indexPatch}(i).level < geometry.maxRefLevel)
                        markRef{indexPatch}(i)=1;
                    end
                    
                    % Get shape functions derivatives
                    [Bu,~,D]=strainGrad(Basis.dgdx{indexPatch}{i},nument,nstress,dim,kgauss,material.C);
                    
                    % Compute strain 
                    strain  = Bu*elemDisp;
                                        
                    % Compute stress                    
                    stress = material.C*strain;
                                            
                    % Compute local stiffness matrix
                    localkUU     = localkUU +gphi.*(Bu'*D).* Basis.volume{indexPatch}{i}(kgauss);
                    
                    % Compute local fintU
                    localfintU    = localfintU + gphi .* Bu' * stress .* Basis.volume{indexPatch}{i}(kgauss);
                    
                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            fintU(elemUDofs)    = fintU(elemUDofs) + localfintU;
            
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
KU = sparse(II,JJ,S,NN,NN);

end




%-------------------------------------------------------------------------%
%  Assemble Phase-field System                                            %
%-------------------------------------------------------------------------%

function [KPf, fintPf] = AssemblePhaseFieldSystem(props,globdat)

Sol        = globdat.state;
oldStepSol = globdat.state0;
PHTelem    = globdat.PHTelem;
MeshInfo   = globdat.mesh;
Basis      = globdat.Basis;

geometry   = props.geom;
material   = props.mat;


% Some useful variables
dim        = geometry.dim;
nstress    = geometry.nstress;
numPatches = geometry.numPatches;
ngaussX    = geometry.ngaussX;
ngaussY    = geometry.ngaussY;
NN         = MeshInfo.sizeBasis;
UEnd       = dim*MeshInfo.sizeBasis;

% Initialise vector for global residual
fintPf        = zeros(NN,1);

% Initialise cells for local stiffness matrices
nelem      = MeshInfo.numElements;
kPhiPhi    = cell(nelem,1);

% Set index counter to zero
elementCounter = 0;
indexCounter   = 0;

% Loop over patches
for indexPatch = 1:numPatches
    
    % Loop over elements
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            
            % Extract element displacement and phase-field
            nument       = size(PHTelem{indexPatch}(i).C,1);
            elemNodes    = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            elemUDofs    = reshape([2*elemNodes-1;2*elemNodes],1,dim*nument);
            elemPhiDofs  = elemNodes+UEnd;
            elemDisp     = Sol(elemUDofs);
            elemPhi      = Sol(elemPhiDofs);
            elemPhi0     = oldStepSol(elemPhiDofs);
            
            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter + nument^2;
            
            % Set local stiffness matrix, residual, displacement to zero
            localkPhiPhi  = zeros(nument,nument);
            localfintPhi  = zeros(nument,1);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    % Get phase-field and compute degradation function
                    phigp   = Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi;
                    phigp0  = Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi0;
                    
                    fac1   = (1-phigp)^material.p;
                    dfac1  = -material.p*(1-phigp)^(material.p-1); 
                    ddfac1 = material.p*(material.p-1)*(1-phigp)^(material.p-2); 
                    
                    fac2   = fac1 + material.a1*phigp + material.a1*material.a2*phigp*phigp + ...
                               material.a1*material.a2*material.a3*phigp*phigp*phigp;
                    dfac2  = dfac1 + material.a1 + 2*material.a1*material.a2*phigp + ...
                               3*material.a1*material.a2*material.a3*phigp*phigp;
                    ddfac2 = ddfac1 + 2*material.a1*material.a2 + ...
                               6*material.a1*material.a2*material.a3*phigp;       
                           
                    % gphi     = fac1/fac2;
                    dgphi    = (dfac1*fac2-fac1*dfac2)/(fac2*fac2);
                    ddgphi   = ((ddfac1*fac2-fac1*ddfac2)*fac2-2*...
                               (dfac1*fac2-fac1*dfac2)*dfac2)/(fac2*fac2*fac2);
                           
                    % alpha  = material.eta*phigp + (1-material.eta)*phigp*phigp                           
                    dalpha   = material.eta + 2*(1-material.eta)*phigp;
                    ddalpha  = 2*(1-material.eta);
                    
                    % Penalisation stuff
                    if phigp - phigp0 > 1e-10
                        delphigp  = 0;
                    else
                        delphigp  = phigp0-phigp;
                    end                    
                    
                    % Get shape functions derivatives
                    [Bu,BPhi,~]=strainGrad(Basis.dgdx{indexPatch}{i},nument,nstress,dim,kgauss,material.C);
                    
                    % Compute gradients (strain and phase-field gradient)
                    strain  = Bu*elemDisp;
                    gradPhi = BPhi*elemPhi;
                                        
                    % Compute stress                    
                    stress = material.C*strain;
                    
                    % Compute new history
                    switch material.Esplit
                        
                        case 'NoSplit'     % Bourdin (2000,2007)
                            Psi = 0.5 * dot(stress,strain);
                            
                        case 'Spectral2D'  % Miehe et. al. (2010)
                            strain_tensor = [strain(1) strain(3)/2; strain(3)/2 strain(2)];  
                            tracep_strain = max(0,strain(1)+strain(2));
                            eigp_strain   = max(0,eig(strain_tensor));
                            Psi           = 0.5 * material.lambda * tracep_strain^2 + ...
                                            material.mu*(eigp_strain(1)^2 + eigp_strain(2)^2); 
                                        
                        case 'Amor'        % Amor et. al. (2009)

                            trueStrain    = Voigt2Tensor(strain,dim,true);
                            trStrainPos   = max(0,trace(trueStrain));
                            volStrain     = trStrainPos/dim;
                            devStrain     = trueStrain-eye(dim)*volStrain;
                            devStrain     = Tensor2Voigt(devStrain,dim);
                            Psi           = 0.5 * material.K * trStrainPos^2 + ...
                                            material.mu * dot(devStrain,devStrain);
                                        
                        case 'Lancioni'    % Lancioni-Royer Carfagni (XXXX)
                            vol_strain    = 0.5*(strain(1)+strain(2));
                            dev_strain    = [strain(1)-vol_strain;
                                             strain(2)-vol_strain;
                                             strain(3)/2];
                            Psi           = material.mu * dot(dev_strain,dev_strain);
                                
                        case 'Rankine'     % Wu (2017)
                            stress_tensor = [stress(1) stress(3)/2; stress(3)/2 stress(2)];
                            eig1_stress   = max(eig(stress_tensor));
                            eig1p_stress  = max(0,eig1_stress);
                            Psi          = 0.5 * eig1p_stress^2/material.E;
                            
                        case 'ModRankine'
                            stress_tensor = [stress(1) stress(3)/2; stress(3) stress(2)/2];
                            eig_stress   = eig(stress_tensor);
                            eigp_stress  = max(0,eig_stress);
                            Psi          = 0.5 * (norm(eigp_stress))^2/material.E;
                        otherwise
                            error('Not implemented!')
                    end

                    % Get initial Psi
                    Psi = max(Psi,material.Psi0);
                        
                    % Compute local stiffness matrix
                    
                    localkPhiPhi = localkPhiPhi + 2*material.gc*material.l0/material.calpha* (BPhi'* BPhi).* Basis.volume{indexPatch}{i}(kgauss);
                    localkPhiPhi = localkPhiPhi + (material.gc/(material.l0*material.calpha)*ddalpha + ddgphi * Psi).*(Basis.shape{indexPatch}{i}(kgauss,:)'* Basis.shape{indexPatch}{i}(kgauss,:)) .* Basis.volume{indexPatch}{i}(kgauss);
                    
                    % Add penalty term
                    
                    if phigp0 > phigp
                        localkPhiPhi = localkPhiPhi + material.pen*material.gc/material.l0.*(Basis.shape{indexPatch}{i}(kgauss,:)'* Basis.shape{indexPatch}{i}(kgauss,:)) .* Basis.volume{indexPatch}{i}(kgauss);
                    end
                    
                    % Compute local residuals
                    
                    localfintPhi  = localfintPhi + 2*material.gc*material.l0/material.calpha* (BPhi'* gradPhi).* Basis.volume{indexPatch}{i}(kgauss);
                    localfintPhi  = localfintPhi + (material.gc/(material.l0*material.calpha)*dalpha + dgphi * Psi - material.pen*material.gc/material.l0*delphigp ) .* Basis.shape{indexPatch}{i}(kgauss,:)' .* Basis.volume{indexPatch}{i}(kgauss);

                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            
            fintPf(elemPhiDofs-UEnd)  =fintPf(elemPhiDofs-UEnd) + localfintPhi;
            
            % Store element matrices in pre-allocated cells
            
            kPhiPhi{elementCounter} = localkPhiPhi;
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
            elmtStiff = [kPhiPhi{elementCounter}];
            II(indexCounter+1:indexCounter+nument^2) = repmat(sctrx,1,nument);
            JJ(indexCounter+1:indexCounter+nument^2) = reshape(repmat(sctrx,nument,1),1,nument^2);
            S(indexCounter+1:indexCounter+nument^2)  = reshape(elmtStiff,1,nument^2);
            indexCounter = indexCounter+nument^2;
        end
    end
end

% Assemble (sparse) stiffness matrix
KPf = sparse(II,JJ,S,NN,NN);

end




%-------------------------------------------------------------------------%
%  Assemble Coupled System                                                %
%-------------------------------------------------------------------------%

function [K, fint, markRef] = AssembleCoupledSystem(props,globdat)

Sol        = globdat.state;
PHTelem    = globdat.PHTelem;
MeshInfo   = globdat.mesh;
Basis      = globdat.Basis;

geometry   = props.geom;
material   = props.mat;


% Some useful variables
dim        = geometry.dim;
nstress    = geometry.nstress;
numPatches = geometry.numPatches;
ngaussX    = geometry.ngaussX;
ngaussY    = geometry.ngaussY;
NN         = globdat.ndofs;
UEnd       = dim*MeshInfo.sizeBasis;

% Initialise cell for mesh refinement marker
markRef    = cell(numPatches,1);

% Initialise vector for global residual
fint        = zeros(NN,1);

% Initialise cells for local stiffness matrices
nelem      = MeshInfo.numElements;
kUU        = cell(nelem,1);
kUPhi      = cell(nelem,1);
kPhiU      = cell(nelem,1);
kPhiPhi    = cell(nelem,1);

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
            elemPhiDofs  = elemNodes+UEnd;
            elemDisp     = Sol(elemUDofs);
            elemPhi      = Sol(elemPhiDofs);
            
            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter + dim^2*nument^2;
            
            % Set local stiffness matrix, residual, displacement to zero
            localkUU      = zeros(dim*nument,dim*nument);
            localkUPhi    = zeros(dim*nument,nument);
            localkPhiPhi  = zeros(nument,nument);
            localfintU    = zeros(dim*nument,1);
            localfintPhi  = zeros(nument,1);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    % Get phase-field and compute degradation function
                    phigp   = Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi;
                    
                    fac1   = (1-phigp)^material.p;
                    dfac1  = -material.p*(1-phigp)^(material.p-1); 
                    ddfac1 = material.p*(material.p-1)*(1-phigp)^(material.p-2); 
                    
                    fac2   = fac1 + material.a1*phigp + material.a1*material.a2*phigp*phigp + ...
                               material.a1*material.a2*material.a3*phigp*phigp*phigp;
                    dfac2  = dfac1 + material.a1 + 2*material.a1*material.a2*phigp + ...
                               3*material.a1*material.a2*material.a3*phigp*phigp;
                    ddfac2 = ddfac1 + 2*material.a1*material.a2 + ...
                               6*material.a1*material.a2*material.a3*phigp;       
                           
                    gphi     = fac1/fac2;
                    dgphi    = (dfac1*fac2-fac1*dfac2)/(fac2*fac2);
                    ddgphi   = ((ddfac1*fac2-fac1*ddfac2)*fac2-2*...
                               (dfac1*fac2-fac1*dfac2)*dfac2)/(fac2*fac2*fac2);
                           
                    % alpha  = material.eta*phigp + (1-material.eta)*phigp*phigp                           
                    dalpha   = material.eta + 2*(1-material.eta)*phigp;
                    ddalpha  = 2*(1-material.eta);
                                       
                    
                    % Set refinement based on phase-field
                    if (phigp > geometry.threshPhi) && (PHTelem{indexPatch}(i).level < geometry.maxRefLevel)
                        markRef{indexPatch}(i)=1;
                    end
                    
                    % Get shape functions derivatives
                    [Bu,BPhi,D]=strainGrad(Basis.dgdx{indexPatch}{i},nument,nstress,dim,kgauss,material.C);
                    
                    % Compute gradients (strain and phase-field gradient)
                    strain  = Bu*elemDisp;
                    gradPhi = BPhi*elemPhi;
                                        
                    % Compute stress                    
                    stress = material.C*strain;
                    
                    % Compute new history
                    switch material.Esplit
                        
                        case 'NoSplit'     % Bourdin (2000,2007)
                            Psi = 0.5 * dot(stress,strain);
                            
                        case 'Spectral2D'  % Miehe et. al. (2010)
                            strain_tensor = [strain(1) strain(3)/2; strain(3)/2 strain(2)];  
                            tracep_strain = max(0,strain(1)+strain(2));
                            eigp_strain   = max(0,eig(strain_tensor));
                            Psi           = 0.5 * material.lambda * tracep_strain^2 + ...
                                            material.mu*(eigp_strain(1)^2 + eigp_strain(2)^2); 
                                        
                        case 'Amor'        % Amor et. al. (2009)

                            trueStrain    = Voigt2Tensor(strain,dim,true);
                            trStrainPos   = max(0,trace(trueStrain));
                            volStrain     = trStrainPos/dim;
                            devStrain     = trueStrain-eye(dim)*volStrain;
                            devStrain     = Tensor2Voigt(devStrain,dim);
                            Psi           = 0.5 * material.K * trStrainPos^2 + ...
                                            material.mu * dot(devStrain,devStrain);
                                        
                        case 'Lancioni'    % Lancioni-Royer Carfagni (XXXX)
                            vol_strain    = 0.5*(strain(1)+strain(2));
                            dev_strain    = [strain(1)-vol_strain;
                                             strain(2)-vol_strain;
                                             strain(3)/2];
                            Psi           = material.mu * dot(dev_strain,dev_strain);
                                
                        case 'Rankine'     % Wu (2017)
                            stress_tensor = [stress(1) stress(3)/2; stress(3)/2 stress(2)];
                            eig1_stress   = max(eig(stress_tensor));
                            eig1p_stress  = max(0,eig1_stress);
                            Psi          = 0.5 * eig1p_stress^2/material.E;
                            
                        case 'ModRankine'
                            stress_tensor = [stress(1) stress(3)/2; stress(3) stress(2)/2];
                            eig_stress   = eig(stress_tensor);
                            eigp_stress  = max(0,eig_stress);
                            Psi          = 0.5 * (norm(eigp_stress))^2/material.E;
                        otherwise
                            error('Not implemented!')
                    end

                    % Get initial Psi
                    Psi = max(Psi,material.Psi0); 
                        
                    % Compute local stiffness matrix
                    localkUU     = localkUU + gphi.*(Bu'*D).* Basis.volume{indexPatch}{i}(kgauss);
                    localkUPhi   = localkUPhi + dgphi .* (Bu' * stress * Basis.shape{indexPatch}{i}(kgauss,:)) .* Basis.volume{indexPatch}{i}(kgauss);
                    localkPhiPhi = localkPhiPhi + 2*material.gc*material.l0/material.calpha* (BPhi'* BPhi).* Basis.volume{indexPatch}{i}(kgauss);
                    localkPhiPhi = localkPhiPhi + (material.gc/(material.l0*material.calpha)*ddalpha + ddgphi * Psi).*(Basis.shape{indexPatch}{i}(kgauss,:)'* Basis.shape{indexPatch}{i}(kgauss,:)) .* Basis.volume{indexPatch}{i}(kgauss);
                    
                    % Compute local residuals
                    localfintU    = localfintU + gphi .* Bu' * stress .* Basis.volume{indexPatch}{i}(kgauss);
                    localfintPhi  = localfintPhi + 2*material.gc*material.l0/material.calpha* (BPhi'* gradPhi).* Basis.volume{indexPatch}{i}(kgauss);
                    localfintPhi  = localfintPhi + (material.gc/(material.l0*material.calpha)*dalpha + dgphi * Psi ) .* Basis.shape{indexPatch}{i}(kgauss,:)' .* Basis.volume{indexPatch}{i}(kgauss);

                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            fint(elemUDofs)    = fint(elemUDofs) + localfintU;
            fint(elemPhiDofs)  = fint(elemPhiDofs) + localfintPhi;
            
            % Store element matrices in pre-allocated cells
            kUU{elementCounter}     = localkUU;
            kUPhi{elementCounter}   = localkUPhi;
            kPhiPhi{elementCounter} = localkPhiPhi;
            kPhiU{elementCounter}   = localkUPhi';
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
            dsctrx    = [reshape([2*sctrx-1;2*sctrx],1,(dim)*nument), UEnd + sctrx];
            elmtStiff = [kUU{elementCounter} kUPhi{elementCounter};kPhiU{elementCounter} kPhiPhi{elementCounter}];
            II(indexCounter+1:indexCounter+(dim+1)^2*nument^2) = repmat(dsctrx,1,(dim+1)*nument);
            JJ(indexCounter+1:indexCounter+(dim+1)^2*nument^2) = reshape(repmat(dsctrx,(dim+1)*nument,1),1,(dim+1)^2*nument^2);
            S(indexCounter+1:indexCounter+(dim+1)^2*nument^2)  = reshape(elmtStiff,1,(dim+1)^2*nument^2);
            indexCounter = indexCounter +(dim+1)^2*nument^2;
        end
    end
end

% Assemble (sparse) stiffness matrix
K = sparse2(II,JJ,S,NN,NN);

end




%-------------------------------------------------------------------------%
% Assemble Fracture Dissipation                                           %
%-------------------------------------------------------------------------%

function [h, g] = AssembleFractureDissipation(props,globdat)

Sol        = globdat.state;
DSol       = globdat.Dstate;
PHTelem    = globdat.PHTelem;
MeshInfo   = globdat.mesh;
Basis      = globdat.Basis;

geometry   = props.geom;
material   = props.mat;

% Some useful variables
dim        = geometry.dim;
numPatches = geometry.numPatches;
nstress    = geometry.nstress;
ngaussX    = geometry.ngaussX;
ngaussY    = geometry.ngaussY;
NN         = (dim+1)*MeshInfo.sizeBasis;
UEnd       = dim*MeshInfo.sizeBasis;

% Initialise h
h = zeros(NN,1);
g = 0;

% Set index counter to zero
elementCounter = 0;
indexCounter   = 0;

% Loop over patches
for indexPatch = 1:numPatches
    
    % Loop over elements
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            
            % Extract element displacement and phase-field
            nument         = size(PHTelem{indexPatch}(i).C,1);
            elemNodes      = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            elemPhiDofs    = elemNodes+UEnd;
            elemPhi        = Sol(elemPhiDofs);
            DelemPhi       = DSol(elemPhiDofs);
            
            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter + dim^2*nument^2;
            
            localh  = zeros(nument,1);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    % Get phase-field and compute degradation function
                    phigp  = Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi;
                    Dphigp = Basis.shape{indexPatch}{i}(kgauss,:) * DelemPhi;
                    
                    % Get shape functions derivatives
                    [~,BPhi,~]=strainGrad(Basis.dgdx{indexPatch}{i},nument,nstress,dim,kgauss,material.C);
                    
                    % Compute phase-field gradient
                    gradPhi    = BPhi*elemPhi;
                    DgradPhi   = BPhi*DelemPhi;
                                       
                    % Compute g
                    g = g + material.gc/material.l0.* phigp*Dphigp.* Basis.volume{indexPatch}{i}(kgauss);
                    g = g + material.gc*material.l0.* dot(gradPhi,DgradPhi).* Basis.volume{indexPatch}{i}(kgauss);
                    
                    % Compute h
                    localh = localh + material.gc/material.l0*(phigp+Dphigp).*Basis.shape{indexPatch}{i}(kgauss,:)'.* Basis.volume{indexPatch}{i}(kgauss);
                    localh = localh + material.gc*material.l0*( BPhi'*(gradPhi+DgradPhi)).* Basis.volume{indexPatch}{i}(kgauss);
                    
                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            h(elemPhiDofs)  = h(elemPhiDofs) + localh;

        end
    end
end

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




%-------------------------------------------------------------------------%
% Assemble Phase-field Mass                                               %
%-------------------------------------------------------------------------%

function M = AssemblePhaseFieldMass(props,globdat)

PHTelem    = globdat.PHTelem;
MeshInfo   = globdat.MeshInfo;
Basis      = globdat.Basis;

geometry   = props.geometry;

% Some useful variables
numPatches = geometry.numPatches;
ngaussX    = geometry.ngaussX;
ngaussY    = geometry.ngaussY;
NN         = MeshInfo.sizeBasis;

% Initialise cells for local stiffness matrices
nelem      = MeshInfo.numElements;
kPhiPhi    = cell(nelem,1);

% Set index counter to zero
elementCounter = 0;
indexCounter   = 0;

% Loop over patches
for indexPatch = 1:numPatches
    
    % Loop over elements
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            
            % Extract element displacement and phase-field
            nument       = size(PHTelem{indexPatch}(i).C,1);            

            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter + nument^2;
            
            % Set local matrix to zero
            localkPhiPhi  = zeros(nument,nument);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                        
                    % Compute local phase-field mass matrix
                    
                    localkPhiPhi = localkPhiPhi + (Basis.shape{indexPatch}{i}(kgauss,:)'* Basis.shape{indexPatch}{i}(kgauss,:)).* Basis.volume{indexPatch}{i}(kgauss);
                    
                end % j Gauss
            end % i Gauss
            
            % Store element matrices in pre-allocated cells            
            kPhiPhi{elementCounter} = localkPhiPhi;
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
            elmtStiff = [kPhiPhi{elementCounter}];
            II(indexCounter+1:indexCounter+nument^2) = repmat(sctrx,1,nument);
            JJ(indexCounter+1:indexCounter+nument^2) = reshape(repmat(sctrx,nument,1),1,nument^2);
            S(indexCounter+1:indexCounter+nument^2)  = reshape(elmtStiff,1,nument^2);
            indexCounter = indexCounter+nument^2;
        end
    end
end

% Assemble (sparse) stiffness matrix
M = sparse(II,JJ,S,NN,NN);

M = M * ones(NN,1);

end



