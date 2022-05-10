function [PhPHTelem,PhcontrolPts,PhmeshInfo,PhmarkRef,Sol,Res,oldSol,oldRes] = RefineTransferNodalDataGlob2Loc(markRef,geometry,Sol,Res,oldSol,oldRes,PhcontrolPts,PhPHTelem,PhmeshInfo)

PhmarkRef = cell(1,geometry.numPatches);                                     

for iPatch = 1:geometry.numPatches
    PhmarkRef(iPatch) = {zeros(length(markRef{iPatch}),1)};
    if (sum(markRef{iPatch})>0)
        
        % Refine and Update the mesh
        disp(['In patch ', num2str(iPatch), ' refining ', num2str(sum(markRef{iPatch})), ' elements.'])
        
        [PhPHTelem{iPatch},PhcontrolPts{iPatch},PhmeshInfo.dimBasis(iPatch),...
         Sol.Phi{iPatch},Sol.U{iPatch},Sol.V{iPatch}, ...
         Res.Phi{iPatch},Res.U{iPatch},Res.V{iPatch}, ...
         oldSol.Phi{iPatch},oldSol.U{iPatch},oldSol.V{iPatch}, ...
         oldRes.Phi{iPatch},oldRes.U{iPatch},oldRes.V{iPatch}, ...
         PhmeshInfo.numElements,PhmarkRef{iPatch}] = ...
            refineElemProjAll(markRef{iPatch},PhPHTelem{iPatch}, ...
            PhcontrolPts{iPatch},geometry,PhmeshInfo.dimBasis(iPatch),...
            Sol.Phi{iPatch},Sol.U{iPatch},Sol.V{iPatch}, ...
            Res.Phi{iPatch},Res.U{iPatch},Res.V{iPatch}, ...
            oldSol.Phi{iPatch},oldSol.U{iPatch},oldSol.V{iPatch}, ...
            oldRes.Phi{iPatch},oldRes.U{iPatch},oldRes.V{iPatch}, ...
            PhmeshInfo.numElements);
    end
end

[PhPHTelem,PhcontrolPts,PhmeshInfo, ...
 Sol.Phi,Sol.U,Sol.V, ...
 Res.Phi,Res.U,Res.V, ...
 oldSol.Phi,oldSol.U,oldSol.V, ...
 oldRes.Phi,oldRes.U,oldRes.V, ...
 PhmarkRef] = checkConformingProjElemAll(PhPHTelem,PhcontrolPts, ...
 PhmeshInfo,geometry, ...
 Sol.Phi,Sol.U,Sol.V, ...
 Res.Phi,Res.U,Res.V, ...
 oldSol.Phi,oldSol.U,oldSol.V, ...
 oldRes.Phi,oldRes.U,oldRes.V, ...
 PhmarkRef);
 
[PhPHTelem,PhmeshInfo] = zipConforming(PhPHTelem,geometry,PhmeshInfo);

SolU   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,Sol.U);
SolV   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,Sol.V);
SolPhi = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,Sol.Phi);

ResU   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,Res.U);
ResV   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,Res.V);
ResPhi = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,Res.Phi);

oldSolU   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,oldSol.U);
oldSolV   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,oldSol.V);
oldSolPhi = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,oldSol.Phi);

oldResU   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,oldRes.U);
oldResV   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,oldRes.V);
oldResPhi = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,oldRes.Phi);

n = PhmeshInfo.sizeBasis;
Sol    = zeros(3*n,1);
Res    = zeros(3*n,1);
oldSol = zeros(3*n,1);
oldRes = zeros(3*n,1);

Sol(1:2:2*n)   = SolU;
Sol(2:2:2*n)   = SolV;
Sol(2*n+1:end) = SolPhi;

Res(1:2:2*n)   = ResU;
Res(2:2:2*n)   = ResV;
Res(2*n+1:end) = ResPhi;

oldSol(1:2:2*n)   = oldSolU;
oldSol(2:2:2*n)   = oldSolV;
oldSol(2*n+1:end) = oldSolPhi;

oldRes(1:2:2*n)   = oldResU;
oldRes(2:2:2*n)   = oldResV;
oldRes(2*n+1:end) = oldResPhi;

end

