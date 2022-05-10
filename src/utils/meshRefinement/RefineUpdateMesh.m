function [PhPHTelem,PhcontrolPts,PhmeshInfo,PhmarkRef] = RefineUpdateMesh(markRef,geometry,solPhiPatch,PhcontrolPts,PhPHTelem,PhmeshInfo)
         
PhmarkRef = cell(1,geometry.numPatches);                                     

for iPatch = 1:geometry.numPatches
    PhmarkRef(iPatch) = {zeros(length(markRef{iPatch}),1)};
    if (sum(markRef{iPatch})>0)
        
        % Refine and Update the mesh
        disp(['In patch ', num2str(iPatch), ' refining ', num2str(sum(markRef{iPatch})), ' elements.'])
        [PhPHTelem{iPatch},PhcontrolPts{iPatch},PhmeshInfo.dimBasis(iPatch),solPhiPatch{iPatch},PhmeshInfo.numElements,PhmarkRef{iPatch}] = ...
            refineElemProj(markRef{iPatch},PhPHTelem{iPatch},PhcontrolPts{iPatch},geometry,PhmeshInfo.dimBasis(iPatch),solPhiPatch{iPatch},PhmeshInfo.numElements);
    end
end

[PhPHTelem,PhcontrolPts,PhmeshInfo,solPhiPatch,PhmarkRef] = checkConformingProjElem(PhPHTelem,PhcontrolPts,PhmeshInfo,geometry,solPhiPatch,PhmarkRef);
[PhPHTelem,PhmeshInfo] = zipConforming(PhPHTelem,geometry,PhmeshInfo);              

% plot1 = subplot(2,2,[1,2]);
% cla(plot1)
% plotMesh(PhPHTelem,PhcontrolPts,geometry,0)
% axis equal

% Phidisp = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,solPhiPatch);

end