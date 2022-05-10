function [PHTelem,meshInfo] = checkConforming_plateWHole(PHTelem,controlPts,geometry,meshInfo)
%checks that the patches are conforming and if needed makes them conforming through mesh refinement
%patchBoundaries format:
% patchA, patchB, edgeA, edgeB
% first patchA should be 1
%edge format: 1-down, 2-right, 3-up, 4-left

p = geometry.p;
q = geometry.q;
patchBoundaries = geometry.patchBoundaries;
numPatches = geometry.numPatches;
numBoundaries = size(patchBoundaries,1);
numberElements = meshInfo.numElements;

keepChecking = 1;
while keepChecking
    keepChecking = 0;
    for patchIndex = 1:numPatches
        for elemIndex = 1:length(PHTelem{patchIndex})
            PHTelem{patchIndex}(elemIndex).nodesGlobal = PHTelem{patchIndex}(elemIndex).nodes;
        end
    end    
    for boundaryIndex = 1:numBoundaries
        % Get the nodes on the boundary edge in patchA and patchB
        patchAList = patchBoundaries{boundaryIndex,1};
        patchB = patchBoundaries{boundaryIndex,2};
        edgeAList = patchBoundaries{boundaryIndex,3};
        edgeBList = patchBoundaries{boundaryIndex,4};
        for indexPatch=1:length(patchAList)
            patchA = patchAList(indexPatch);
            edgeA = edgeAList(indexPatch);
            edgeB = edgeBList(indexPatch);
            quadListA = meshInfo.quadList{patchA};
            quadListB = meshInfo.quadList{patchB};
        
            [elementsA] = sortEdgeElem( PHTelem{patchA}, edgeA);
            [elementsB] = sortEdgeElem( PHTelem{patchB}, edgeB);
        
            [quadRefA, quadRefB] = makeConforming(PHTelem{patchA}, PHTelem{patchB}, elementsA{1}, elementsB{1}, edgeA, edgeB, quadListA, quadListB);
            indexQuadA = find(quadRefA > 0);
            indexQuadB = find(quadRefB > 0);
        
            if ~isempty(indexQuadA)
                keepChecking = 1;
                numNewPatches = length(indexQuadA);
                disp(['In patch ', num2str(patchA), ' refining ',num2str(numNewPatches), ' quadruplets to keep conformity with patch ', num2str(patchB)])
                [meshInfo.quadList{patchA},PHTelem{patchA},controlPts{patchA},meshInfo.dimBasis(patchA)] = refineMeshIso(quadRefA,quadListA,PHTelem{patchA},controlPts{patchA},p,q,dimBasis(patchA));
                numberElements = numberElements + (size(meshInfo.quadList{patchA},1) - size(quadListA,1))*4;
            end
        
            if ~isempty(indexQuadB)
                keepChecking = 1;
                numNewPatches = length(indexQuadB);
                disp(['In patch ', num2str(patchB), ' refining ',num2str(numNewPatches), ' quadruplets to keep conformity with patch ', num2str(patchA)])
                [meshInfo.quadList{patchB},PHTelem{patchB},controlPts{patchB},meshInfo.dimBasis(patchB)] = refineMeshIso(quadRefB,quadListB,PHTelem{patchB},controlPts{patchB},p,q,dimBasis(patchB));
                numberElements = numberElements + (size(meshInfo.quadList{patchB},1) - size(quadListB,1))*4;
            end
        end
    end
end

meshInfo.numElements = numberElements;
end