function [PHTelem,controlPts,dimBasis,quadList,markRef] = checkConformingMark_plateWHole(PHTelem,controlPts,dimBasis,patchBoundaries,p,q,quadList,markRef)
%checks that the patches are conforming and if needed makes them conforming through mesh refinement
%patchBoundaries format:
% patchA, patchB, edgeA, edgeB
% first patchA should be 1
%edge format: 1-down, 2-right, 3-up, 4-left

numBoundaries = size(patchBoundaries,1);
numPatches = length(PHTelem);

%create/set nodesGlobal entries in all patches to be equal to local nodes
%entries
keepChecking = 1;
while keepChecking
    keepChecking = 0;
    for patchIndex = 1:numPatches
        for elemIndex = 1:length(PHTelem{patchIndex})
            PHTelem{patchIndex}(elemIndex).nodesGlobal = PHTelem{patchIndex}(elemIndex).nodes;
        end
    end
    
    for boundaryIndex = 1:numBoundaries
        %get the nodes on the boundary edge in patchA and patchB
        patchAList = patchBoundaries{boundaryIndex,1};
        patchB = patchBoundaries{boundaryIndex,2};
        edgeAList = patchBoundaries{boundaryIndex,3};
        edgeBList = patchBoundaries{boundaryIndex,4};
        
        for indexPatch=1:length(patchAList)
            patchA = patchAList(indexPatch);
            edgeA = edgeAList(indexPatch);
            edgeB = edgeBList(indexPatch);
            quadListA = quadList{patchA};
            quadListB = quadList{patchB};
            
            [elementsA] = sortEdgeElem(PHTelem{patchA}, edgeA);
            [elementsB] = sortEdgeElem(PHTelem{patchB}, edgeB);
            
            [quadRefA, quadRefB] = makeConforming(PHTelem{patchA}, PHTelem{patchB}, elementsA{1}, elementsB{1}, edgeA, edgeB, quadListA, quadListB);
            indexQuadA = find(quadRefA > 0);
            indexQuadB = find(quadRefB > 0);
            
            if ~isempty(indexQuadA)
                keepChecking = 1;
                numNewPatches = length(indexQuadA);
                disp(['In patch ', num2str(patchA), ' refining ',num2str(numNewPatches), ' quadruplets to keep conformity with patch ', num2str(patchB)])
                [quadList{patchA}, PHTelem{patchA}, controlPts{patchA}, dimBasis(patchA),markRefTemp] = refineMeshIsoMark(quadRefA, quadListA, PHTelem{patchA}, controlPts{patchA}, p, q, dimBasis(patchA));
                indexRef = find(markRefTemp > 0);
                if indexRef
                    markRef{patchA}(indexRef) = 1;
                end
            end
            
            if ~isempty(indexQuadB)
                keepChecking = 1;
                numNewPatches = length(indexQuadB);
                disp(['In patch ', num2str(patchB), ' refining ',num2str(numNewPatches), ' quadruplets to keep conformity with patch ', num2str(patchA)])
                [quadList{patchB}, PHTelem{patchB}, controlPts{patchB}, dimBasis(patchB),markRefTemp] = refineMeshIsoMark(quadRefB, quadListB, PHTelem{patchB}, controlPts{patchB}, p, q, dimBasis(patchB));
                indexRef = find(markRefTemp > 0);
                if indexRef
                    markRef{patchB}(indexRef) = 1;
                end
            end
        end
    end
end