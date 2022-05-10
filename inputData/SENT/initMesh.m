function [PHTelem,controlPts,meshInfo] = initMesh(geometry)

meshInfo = struct;
L = geometry.L;
W = geometry.W;
p = geometry.p;
q = geometry.q;
numberElemU = geometry.numElemU;
numberElemV = geometry.numElemV;
numPatches = geometry.numPatches;

% Initialize the PHT geometry on coarse mesh
controlPts = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
dimBasis = zeros(1, numPatches);
quadList = cell(numPatches,1);

% Divide the patches along the x direction
xVertices = linspace(0,L,3);
yVertices = linspace(0,W,3);

patchCounter = 0;
patchIndexSet = [1,2,4,3];
numberElements = 0;

for patchIndexY = 1:(numPatches-2)
    for patchIndexX = 1:(numPatches-2)
        % Set the dimensions of the patch
        patchCounter = patchCounter + 1;
        patchIndex = patchIndexSet(patchCounter);
        patchMinX = xVertices(patchIndexX);
        patchMaxX = xVertices(patchIndexX+1);
        patchMinY = yVertices(patchIndexY);
        patchMaxY = yVertices(patchIndexY+1);
        
        % Initialize geometry on coarsest mesh
        coefs(1:3,1,1) = [patchMinX; patchMinY;0];
        coefs(1:3,1,2) = [patchMinX; patchMaxY;0];
        coefs(1:3,2,1) = [patchMaxX; patchMinY;0];
        coefs(1:3,2,2) = [patchMaxX; patchMaxY;0];
        coefs(4,1,1) = 1;
        coefs(4,1,2) = 1;
        coefs(4,2,1) = 1;
        coefs(4,2,2) = 1;
        
        knotU = [0 0 1 1];
        knotV = [0 0 1 1];
        
        nurbs = nrbmak(coefs,{knotU,knotV});
        p_init = nurbs.order(1)-1;
        q_init = nurbs.order(2)-1;
        
        % Refine into numberElemU by numberElemV knotspans
        knotU = linspace(0,1,numberElemU+1);
        knotV = linspace(0,1,numberElemV+1);
        
        numberElementsU = length(unique(knotU))-1;
        numberElementsV = length(unique(knotV))-1;
        
        % Increase polynomial order
        nurbs = nrbdegelev(nurbs,[p-p_init,q-q_init]);
        nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});
        % Repeat the knots to get C1 continuity
        nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});

        [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = genControlPtsNoRep(nurbs,p,q,numberElementsU,numberElementsV);
        
        % Refine the area between (0.45...0.55) in the y-direction
        [quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,0.45+0.01,0.55-0.01,0.45+0.01,0.55-0.01);
        [quadList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef,quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},p,q,dimBasis(patchIndex));
        
        % Refine the area between (0.495...0.505) in the y-direction
        [quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,0.475+0.001,0.525-0.001,0.475+0.001,0.525-0.001);
        [quadList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef,quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},p,q,dimBasis(patchIndex));
        
        % Refine the area between (0.4995...0.5005) in the y-direction
        [quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,0.4875+0.0001,0.5125-0.0001,0.4875+0.0001,0.5125-0.0001);
        [quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef, quadList{patchIndex}, PHTelem{patchIndex}, controlPts{patchIndex},p,q,dimBasis(patchIndex));
       
        numberElements = numberElements + 4*size(quadList{patchIndex},1);
    end
end

meshInfo.dimBasis = dimBasis;
meshInfo.quadList = quadList;
meshInfo.numElements = numberElements;

end