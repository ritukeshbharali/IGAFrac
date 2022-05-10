function [PhPHTelem,PhcontrolPts,PhmeshInfo,PhmarkRef,field] = RefineTransferNodalLoc2Glob(PhPHTelem,PhcontrolPts,PhmeshInfo,markRef,geometry,field)

PhmarkRef = cell(1,geometry.numPatches);                                     

for iPatch = 1:geometry.numPatches
    PhmarkRef(iPatch) = {zeros(length(markRef{iPatch}),1)};
    if (sum(markRef{iPatch})>0)
        
        % Refine and Update the mesh
        disp(['In patch ', num2str(iPatch), ' refining ', num2str(sum(markRef{iPatch})), ' elements.'])
        
        [PhPHTelem{iPatch},PhcontrolPts{iPatch},PhmeshInfo.dimBasis(iPatch),...
         field, ...
         PhmeshInfo.numElements,PhmarkRef{iPatch}] = ...
            refineElemProjField(markRef{iPatch},PhPHTelem{iPatch}, ...
            PhcontrolPts{iPatch},geometry,PhmeshInfo.dimBasis(iPatch),...
            geometry.dim,field, ...
            PhmeshInfo.numElements);
    end
end

[PhPHTelem,PhcontrolPts,PhmeshInfo, ...
 field, ...
 PhmarkRef] = checkConformingProjElemAll(PhPHTelem,PhcontrolPts, ...
 PhmeshInfo,geometry, ...
 field, ...
 PhmarkRef);
 
[PhPHTelem,PhmeshInfo] = zipConforming(PhPHTelem,geometry,PhmeshInfo);

switch geometry.dim
    
    case 2
        
        fieldX   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,field.X);
        fieldY   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,field.Y);
        fieldPhi = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,field.Phi);
        
        n = PhmeshInfo.sizeBasis;
        field  = zeros(3*n,1);

        field(1:2:2*n)   = fieldX;
        field(2:2:2*n)   = fieldY;
        field(2*n+1:end) = fieldPhi;
        
    case 3
        
        fieldX   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,field.X);
        fieldY   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,field.Y);
        fieldZ   = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,field.Z);
        fieldPhi = transferFieldLoc2Glob(PhPHTelem,PhmeshInfo,field.Phi);
        
        n = PhmeshInfo.sizeBasis;
        field  = zeros(4*n,1);

        field(1:2:3*n)   = fieldX;
        field(2:2:3*n)   = fieldY;
        field(3:2:3*n)   = fieldZ;
        field(3*n+1:4*n) = fieldPhi;
        
end

end

