function [fieldData] = transferFieldLoc2Glob(PHTelem,PhmeshInfo,fieldDataPatch)
% Transfer the fieldData indices from the global array to the local patch
% structure

numPatches = length(PHTelem);
fieldData = zeros(PhmeshInfo.sizeBasis,size(fieldDataPatch{1},2));

for iPatch = 1:numPatches
    for iElem = 1:length(PHTelem{iPatch})
        if isempty(PHTelem{iPatch}(iElem).children)
            globIndex = PHTelem{iPatch}(iElem).nodesGlobal;
            localIndex = PHTelem{iPatch}(iElem).nodes;
            fieldData(globIndex,:) = fieldDataPatch{iPatch}(localIndex,:);
        end
    end
end