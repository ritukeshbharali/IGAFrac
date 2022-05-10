function [PHTelem,meshInfo] = zipConforming(PHTelem,geometry,meshInfo)
% Connects two conforming patches by changing the nodesGlobal entry
% patchBoundaries format:
% patchA, patchB, edgeA, edgeB, flagRepeat
% patchA should be 1
% flagRepeat = 0 --> shift node indices, flagRepeat = 1 -->
% edge format: 1-down, 2-right, 3-up, 4-left

p = geometry.p;
q = geometry.q;
patchBoundaries = geometry.patchBoundaries;
numPatches = geometry.numPatches;
numBoundaries = size(patchBoundaries,1);

%create/set nodesGlobal entries in all patches to be equal to local nodes
%entries
for patchIndex = 1:numPatches
    for elemIndex = 1:length(PHTelem{patchIndex})
        PHTelem{patchIndex}(elemIndex).nodesGlobal = PHTelem{patchIndex}(elemIndex).nodes;
    end
end

curShift = meshInfo.dimBasis(1);
overlapCounter = 0;
patchesSeen = [];

for boundaryIndex = 1:numBoundaries
    
    % Get the nodes on the boundary edge in patchA and patchB
    patchAList = patchBoundaries{boundaryIndex,1};
    patchB = patchBoundaries{boundaryIndex,2};
    edgeAList = patchBoundaries{boundaryIndex,3};
    edgeBList = patchBoundaries{boundaryIndex,4};
    nodesPattern = zeros(1,meshInfo.dimBasis(patchB));
    nodesA = [];
    nodesB = [];
    patchesSeen = union(patchesSeen, patchB);
    
    for indexPatch=1:length(patchAList)
        patchA = patchAList(indexPatch);
        edgeA = edgeAList(indexPatch);
        edgeB = edgeBList(indexPatch);
        patchesSeen = union(patchesSeen, patchA);
        
        nodesAcell = sortEdgeNodesElem( PHTelem{patchA}, edgeA, p, q );
        nodesBcell = sortEdgeNodesElem( PHTelem{patchB}, edgeB, p, q );
        
        nodesA = [nodesA, nodesAcell{1}];
        nodesB = [nodesB, nodesBcell{1}];
        
        if length(nodesA)~=length(nodesB)
            error('Non-conforming patches encountered. Aborting...')
        end
    end
    [nodesB,sI] = sort(nodesB);
    nodesA = nodesA(sI);

    curBdryNode = 0;
    for nodeIndex=1:length(nodesA)
        %shift the basis functions indices in nodesPattern
        prevBdryNode = curBdryNode;
        curBdryNode = nodesB(nodeIndex);
        nodesPattern(prevBdryNode+1:curBdryNode-1) = ((prevBdryNode+1):(curBdryNode-1)) + curShift;
        nodesPattern(curBdryNode) = nodesA(nodeIndex);
        if prevBdryNode<curBdryNode
            curShift = curShift - 1;
        end
    end
    %shift the indices after the last boundary node
    nodesPattern(curBdryNode+1:end) = ((curBdryNode+1):meshInfo.dimBasis(patchB)) + curShift;
    %update the nodesGlobal in patchB according to nodesPattern
    for elemIndex = 1:length(PHTelem{patchB})
        PHTelem{patchB}(elemIndex).nodesGlobal = nodesPattern(PHTelem{patchB}(elemIndex).nodes);
    end    
    overlapCounter = overlapCounter + length(unique(nodesA));
    
    curShift = sum(meshInfo.dimBasis(patchesSeen))-overlapCounter;    
    
end

meshInfo.sizeBasis = sum(meshInfo.dimBasis) - overlapCounter;

end



