function [PhPHTelem,EPHTelem,PhcontrolPts,EcontrolPts,PhmeshInfo,EmeshInfo,Phidisp,PhmarkRef,EmarkRef] = ...
    refineMesh(PhPHTelem,PhmeshInfo,Phidisp,geometry,PhcontrolPts,EPHTelem,EmeshInfo,EcontrolPts,EoctupleRef,PhmarkRef)

[solPhiPatch] = transferFieldGlob2Loc(PhPHTelem,PhmeshInfo,Phidisp);
PhdimBasis = PhmeshInfo.dimBasis;
PhNumElements = PhmeshInfo.numElements;
parfor iPatch = 1:geometry.numPatches
    if (sum(PhmarkRef{iPatch})>0)
        % Refine and Update the mesh
        disp(['In patch ', num2str(iPatch), ' refining ', num2str(sum(PhmarkRef{iPatch})), ' elements...'])
        [PhPHTelem{iPatch},PhcontrolPts{iPatch},PhdimBasis(iPatch),solPhiPatch{iPatch},PhNumElements(iPatch),PhmarkRef{iPatch}] = ...
            refineElemProjGradedIso3D(PhmarkRef{iPatch},PhPHTelem{iPatch},PhcontrolPts{iPatch},geometry,PhdimBasis(iPatch),...
            solPhiPatch{iPatch},PhNumElements(iPatch));
    end
end

[PhPHTelem,PhcontrolPts,PhmeshInfo.dimBasis,solPhiPatch,PhmeshInfo.numElements,PhmarkRef] = ...
    checkConformingProjElemMark3D(PhPHTelem,PhcontrolPts,PhdimBasis,geometry,solPhiPatch,PhNumElements,PhmarkRef);
[PhPHTelem,PhmeshInfo] = zipConforming3D(PhPHTelem,PhmeshInfo,geometry);
Phidisp = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,solPhiPatch);

EmarkRef = cell(1,geometry.numPatches);
EdimBasis = EmeshInfo.dimBasis;
EOctupleList = EmeshInfo.octupleList;
if sum([EoctupleRef{:}])
    for iPatch = 1:geometry.numPatches
        EmarkRef(iPatch) = {zeros(length(EPHTelem{iPatch}),1)};
        if sum(EoctupleRef{iPatch})>0
            
            [EOctupleList{iPatch},EPHTelem{iPatch},EcontrolPts{iPatch},EdimBasis(iPatch),EmarkRef{iPatch}] = ...
                refineMeshGraded3D(EoctupleRef{iPatch},EOctupleList{iPatch},EPHTelem{iPatch},EcontrolPts{iPatch},geometry,EdimBasis(iPatch));
        end
    end
    clear EoctupleRef
    [EPHTelem,EcontrolPts,EmeshInfo.dimBasis,EmeshInfo.octupleList,EmarkRef] = checkConformingMark3D(EPHTelem,EcontrolPts,EdimBasis,geometry,EOctupleList,EmarkRef);
    [EPHTelem,EmeshInfo] = zipConforming3D(EPHTelem,EmeshInfo,geometry);
end
end