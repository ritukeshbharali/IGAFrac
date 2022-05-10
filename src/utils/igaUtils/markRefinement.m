function [PhmarkRef,EmarkRef] = markRefinement(PhmarkRef,EmarkRef,ERefelem,PhRefelem)
% Obtains the elements to be refined on the phase field mesh and the
% elastic field mesh

numPatches =  length(PhmarkRef);
tempEmark = cell(numPatches,1);
tempPhmark = cell(numPatches,1);

for indexPatch = 1:numPatches
    tempEmark(indexPatch) = {zeros(length(EmarkRef{indexPatch}),1)};
    tempPhmark(indexPatch) = {zeros(length(PhmarkRef{indexPatch}),1)};
    ELength = length(ERefelem{indexPatch});% checking if conforming algorithm works, the elements take the children of the neighbouring elements which are refined in that particular iteration
    for i = 1:length(PhmarkRef{indexPatch})
        if (PhmarkRef{indexPatch}(i))&& (i <= ELength)
            tempPhmark{indexPatch}(i) = 1;
            ElasElmtNum = unique(ERefelem{indexPatch}{i}); % Locating the element coresponding to element 'i' in the Phase Mesh to the Elastic mesh
            tempEmark{indexPatch}(ElasElmtNum) = 1;
        end
    end
    PhLength = length(PhRefelem{indexPatch});
    for i = 1:length(EmarkRef{indexPatch})
        if (EmarkRef{indexPatch}(i)) && (i <= PhLength)
            tempEmark{indexPatch}(i) = 1;
            PhaseElmtNum = unique(PhRefelem{indexPatch}{i}); % Locating the element coresponding to element 'i' in the Phase Mesh to the Elastic mesh
            tempPhmark{indexPatch}(PhaseElmtNum) = 1;
        end
    end
end

PhmarkRef = tempPhmark;
EmarkRef = tempEmark;

end