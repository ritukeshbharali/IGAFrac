function [PHTelem,controlPts,dimBasis,numberElements,markRef] = checkConformingElemIsoMark(PHTelem,controlPts,dimBasis,patchBoundaries,p,q,numberElements,markRef)
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
            [elementsA] = sortEdgeElem( PHTelem{patchA}, edgeA);
            [elementsB] = sortEdgeElem( PHTelem{patchB}, edgeB);
            
            [elemRefA, elemRefB] = makeConformingElem(PHTelem{patchA}, PHTelem{patchB}, elementsA{1}, elementsB{1}, edgeA, edgeB);
            indexElemA = find(elemRefA > 0);
            indexElemB = find(elemRefB > 0);
            
            if ~isempty(indexElemA)
                keepChecking = 1;
                numNewElems = length(indexElemA);
                disp(['In patch ', num2str(patchA), ' refining ',num2str(numNewElems), ' elements to keep conformity with patch ', num2str(patchB)])
                [PHTelem{patchA},controlPts{patchA},dimBasis(patchA),numberElements,markRefTemp] = refineElemGradedIso(elemRefA,PHTelem{patchA},controlPts{patchA},p,q,dimBasis(patchA),numberElements);
                indexRef = find(markRefTemp > 0);
                if indexRef
                    markRef{patchA}(indexRef) = 1;
                end
            end
            
            if ~isempty(indexElemB)
                keepChecking = 1;
                numNewElems = length(indexElemB);
                disp(['In patch ', num2str(patchB), ' refining ',num2str(numNewElems), ' elements to keep conformity with patch ', num2str(patchA)])
                [PHTelem{patchB},controlPts{patchB},dimBasis(patchB),numberElements,markRefTemp] = refineElemGradedIso(elemRefB,PHTelem{patchB},controlPts{patchB},p,q,dimBasis(patchB),numberElements);
                indexRef = find(markRefTemp > 0);
                if indexRef
                    markRef{patchB}(indexRef) = 1;
                end
            end
        end
    end
end


