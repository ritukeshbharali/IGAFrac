function [Data11,Data22,Data33,Data44] = TransferNodalDataGlob2Loc(Data1,Data2,Data3,Data4,PHTelem,MeshInfo,dim)

if dim ~= 2
    error('Only 2D supported')
end

% Transfer: Global to Local
Data11.U    = transferFieldGlob2Loc(PHTelem,MeshInfo,Data1(1:2:2*MeshInfo.sizeBasis));
Data11.V    = transferFieldGlob2Loc(PHTelem,MeshInfo,Data1(2:2:2*MeshInfo.sizeBasis));
Data11.Phi  = transferFieldGlob2Loc(PHTelem,MeshInfo,Data1(2*MeshInfo.sizeBasis+1:3*MeshInfo.sizeBasis));
   
Data22.U    = transferFieldGlob2Loc(PHTelem,MeshInfo,Data2(1:2:2*MeshInfo.sizeBasis));
Data22.V    = transferFieldGlob2Loc(PHTelem,MeshInfo,Data2(2:2:2*MeshInfo.sizeBasis));
Data22.Phi  = transferFieldGlob2Loc(PHTelem,MeshInfo,Data2(2*MeshInfo.sizeBasis+1:3*MeshInfo.sizeBasis));
    
Data33.U    = transferFieldGlob2Loc(PHTelem,MeshInfo,Data3(1:2:2*MeshInfo.sizeBasis));
Data33.V    = transferFieldGlob2Loc(PHTelem,MeshInfo,Data3(2:2:2*MeshInfo.sizeBasis));
Data33.Phi  = transferFieldGlob2Loc(PHTelem,MeshInfo,Data3(2*MeshInfo.sizeBasis+1:3*MeshInfo.sizeBasis));
    
Data44.U    = transferFieldGlob2Loc(PHTelem,MeshInfo,Data4(1:2:2*MeshInfo.sizeBasis));
Data44.V    = transferFieldGlob2Loc(PHTelem,MeshInfo,Data4(2:2:2*MeshInfo.sizeBasis));
Data44.Phi  = transferFieldGlob2Loc(PHTelem,MeshInfo,Data4(2*MeshInfo.sizeBasis+1:3*MeshInfo.sizeBasis));
    
end

