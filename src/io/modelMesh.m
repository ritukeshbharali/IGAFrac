function [PHTelem,controlPts,meshInfo] = modelMesh(geometry)
% Generates the initial coarse mesh
% Non-uniform refinement for cubic PHT splines

[PHTelem,controlPts,meshInfo] = initMesh(geometry);
[PHTelem,meshInfo] = checkConforming(PHTelem,geometry,meshInfo);
[PHTelem,meshInfo] = zipConforming(PHTelem,geometry,meshInfo);
    
end