function [PHTelem,controlPts,meshInfo,fieldDataPatch,markRef] = checkConformingProjElem(PHTelem,controlPts,meshInfo,geometry,fieldDataPatch,markRef)
% Checks that the patches are conforming and if needed makes them conforming through mesh refinement
% patchBoundaries format: patchA, patchB, edgeA, edgeB
% patchA should be 1
% edge format: 1-down, 2-right, 3-up, 4-left
% Transfer field data at patch level

patchBoundaries = geometry.patchBoundaries;
numBoundaries = size(patchBoundaries,1);
numPatches = geometry.numPatches;

keepChecking = 1;
while keepChecking
    keepChecking = 0;
    % Create/set nodesGlobal entries in all patches to be equal to local nodes
    % entries
    for patchIndex = 1:numPatches
        for elemIndex = 1:length(PHTelem{patchIndex})
            PHTelem{patchIndex}(elemIndex).nodesGlobal = PHTelem{patchIndex}(elemIndex).nodes;
        end
    end
    
    %loop over the boundaries and match patches if needed
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
            
            [elementsA] = sortEdgeElem( PHTelem{patchA}, edgeA);
            [elementsB] = sortEdgeElem( PHTelem{patchB}, edgeB);
            
            [elemRefA, elemRefB] = makeConformingElem(PHTelem{patchA}, PHTelem{patchB}, elementsA{1}, elementsB{1}, edgeA, edgeB);
            indexElemA = find(elemRefA > 0);
            indexElemB = find(elemRefB > 0);
            
            if ~isempty(indexElemA)
                keepChecking = 1;
                numNewElems = length(indexElemA);
                disp(['In patch ', num2str(patchA), ' refining ',num2str(numNewElems), ' elements to keep conformity with patch ', num2str(patchB)])
                [PHTelem{patchA},controlPts{patchA},meshInfo.dimBasis(patchA),fieldDataPatch{patchA},meshInfo.numberElements,markRefTemp] = ...
                    refineElemProj(elemRefA, PHTelem{patchA},controlPts{patchA},geometry,meshInfo.dimBasis(patchA),fieldDataPatch{patchA},meshInfo.numElements);
                indexRef = find(markRefTemp > 0);
                if indexRef
                    markRef{patchA}(indexRef) = 1;
                end
            end
            
            if ~isempty(indexElemB)
                keepChecking = 1;
                numNewElems = length(indexElemB);
                disp(['In patch ', num2str(patchB), ' refining ',num2str(numNewElems), ' elements to keep conformity with patch ', num2str(patchA)])
                [PHTelem{patchB},controlPts{patchB},meshInfo.dimBasis(patchB),fieldDataPatch{patchB},meshInfo.numberElements,markRefTemp] = ...
                    refineElemProj(elemRefB, PHTelem{patchB},controlPts{patchB},geometry,meshInfo.dimBasis(patchB),fieldDataPatch{patchB},meshInfo.numElements);
                indexRef = find(markRefTemp > 0);
                if indexRef
                    markRef{patchB}(indexRef) = 1;
                end
            end
        end
    end
end
end