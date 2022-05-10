function [patchList,PHTelem,controlPts,dimBasis,markRef] = refineMeshIsoMark(patchRef,patchList,PHTelem,controlPts,p,q,dimBasis)
% Refines the patches marked by patchRef
% Also marks the element that are refined
numElem = length(PHTelem);
indexPatch = find(patchRef > 0);
numNewPatches = length(indexPatch);
newPatchList = zeros(size(patchList,1)+numNewPatches,4);
markRef = zeros(length(PHTelem),1);

%sort patchList by level of elements to prevent coarse->fine level errors
levelList = zeros(1,size(patchList,1));
for i=1:size(patchList,1)
    levelList(i) = PHTelem(patchList(i,1)).level;
end
    
[~,sortedIndex] = sort(levelList);
patchList = patchList(sortedIndex, :);
patchRef = patchRef(sortedIndex);
tempPHTelem = PHTelem;
patchCounter = 0;
for i=1:length(patchRef)
    if patchRef(i)
        curPatchElem = patchList(i,:);
        markRef(curPatchElem) = 1;
        for iCount = 1:4
            curElem = curPatchElem(iCount);
            north = tempPHTelem(curElem).neighbor_up;
            south = tempPHTelem(curElem).neighbor_down;
            east = tempPHTelem(curElem).neighbor_right;
            west = tempPHTelem(curElem).neighbor_left;
            if ~isempty(north)
            northeast = tempPHTelem(max(north)).neighbor_right;
            northwest = tempPHTelem(min(north)).neighbor_left;
            else
                northeast = [];
                northwest = [];
            end
            if ~isempty(south)
            southeast = tempPHTelem(max(south)).neighbor_right;
            southwest = tempPHTelem(min(south)).neighbor_left;
            else
                southeast = [];
                southwest = [];
            end
            neighbor = [east,west,north,south,northeast,northwest,southeast,southwest];  
            markRef(neighbor) = 1;
        end
        [PHTelem,controlPts,dimBasis] = crossInsert(PHTelem,controlPts,curPatchElem,dimBasis,p,q);
        newPatchList(patchCounter+1,:) = numElem+1:numElem+4;
        newPatchList(patchCounter+2,:) = numElem+5:numElem+8;
        newPatchList(patchCounter+3,:) = numElem+9:numElem+12;
        newPatchList(patchCounter+4,:) = numElem+13:numElem+16;
        numElem = numElem + 16;
        patchCounter = patchCounter + 4;
    else
        newPatchList(patchCounter+1,:) = patchList(i,:);
        patchCounter = patchCounter + 1;
    end                    
end
 patchList = newPatchList;
end
