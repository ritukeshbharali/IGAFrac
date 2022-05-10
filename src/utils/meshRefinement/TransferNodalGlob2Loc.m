function field = TransferNodalGlob2Loc(field,PHTelem,MeshInfo,dim)

if dim ~= 2
    error('Only 2D supported')
end

switch dim
    
    case 2
        
        % Transfer: Global to Local
        field      = struct;
        field.X    = transferFieldGlob2Loc(PHTelem,MeshInfo,field(1:2:2*MeshInfo.sizeBasis));
        field.Y    = transferFieldGlob2Loc(PHTelem,MeshInfo,field(2:2:2*MeshInfo.sizeBasis));
        field.Phi  = transferFieldGlob2Loc(PHTelem,MeshInfo,field(2*MeshInfo.sizeBasis+1:3*MeshInfo.sizeBasis));
        
    case 3
        
        % Transfer: Global to Local
        field.X    = transferFieldGlob2Loc(PHTelem,MeshInfo,field(1:2:3*MeshInfo.sizeBasis));
        field.Y    = transferFieldGlob2Loc(PHTelem,MeshInfo,field(2:2:3*MeshInfo.sizeBasis));
        field.Z    = transferFieldGlob2Loc(PHTelem,MeshInfo,field(3:2:3*MeshInfo.sizeBasis));
        field.Phi  = transferFieldGlob2Loc(PHTelem,MeshInfo,field(3*MeshInfo.sizeBasis+1:4*MeshInfo.sizeBasis));
        
    otherwise
        error('Check problem dimensions!')
end

   
end

