function [PhPHTelem,PhcontrolPts,PhmeshInfo,Sol,PhmarkRef] = refinePhaseAll(markRef,geometry,solPhiPatch,solUPatch,solVPatch,PhcontrolPts,PhPHTelem,PhmeshInfo)
         
PhmarkRef = cell(1,geometry.numPatches);                                     

for iPatch = 1:geometry.numPatches
    PhmarkRef(iPatch) = {zeros(length(markRef{iPatch}),1)};
    if (sum(markRef{iPatch})>0)
        
        % Refine and Update the mesh
        disp(['In patch ', num2str(iPatch), ' refining ', num2str(sum(markRef{iPatch})), ' elements.'])
        [PhPHTelem{iPatch},PhcontrolPts{iPatch},PhmeshInfo.dimBasis(iPatch),solPhiPatch{iPatch},solUPatch{iPatch},solVPatch{iPatch},PhmeshInfo.numElements,PhmarkRef{iPatch}] = ...
            refineElemProjAll(markRef{iPatch},PhPHTelem{iPatch},PhcontrolPts{iPatch},geometry,PhmeshInfo.dimBasis(iPatch),solPhiPatch{iPatch},solUPatch{iPatch},solVPatch{iPatch},PhmeshInfo.numElements);
    end
end

[PhPHTelem,PhcontrolPts,PhmeshInfo,solPhiPatch,solUPatch,solVPatch,PhmarkRef] = checkConformingProjElemAll(PhPHTelem,PhcontrolPts,PhmeshInfo,geometry,solPhiPatch,solUPatch,solVPatch,PhmarkRef);
[PhPHTelem,PhmeshInfo] = zipConforming(PhPHTelem,geometry,PhmeshInfo);

% At this point, the mesh is refined. The number of elements per patch
% increases. For solPhiPatch, the crossInsertProjIso has been carried out.
% But solUPatch and solVPatch still has the old element numbers
% (crossInsertProjIso not done!). So trasferFieldLoc2Glob fails (array
% index exceeded error).

% plot1 = subplot(2,2,[1,2]);
% cla(plot1)
% plotMesh(PhPHTelem,PhcontrolPts,geometry,0)
% axis equal

Phi = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,solPhiPatch);
U   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,solUPatch);
V   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,solVPatch);

n = PhmeshInfo.sizeBasis;
Sol = zeros(3*n,1);
Sol(1:2:2*n)   = U;
Sol(2:2:2*n)   = V;
Sol(2*n+1:end) = Phi;

end