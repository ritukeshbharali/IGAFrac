function [elemRefA,elemRefB] = makeConformingElem(PHTelemA,PHTelemB,elementsA,elementsB,edgeA,edgeB)
% Checks that patchA and patchB are conforming (the elements match up) and
% refines if necessary
% Note: assumes only 1 level of refinement is needed
% operates at element rather than patch level


%initialize elemRefA and elemRefB with zeros (no refinements)
elemRefA = zeros(1, length(PHTelemA));
elemRefB = zeros(1, length(PHTelemB));

%make list of coordinates of element vertices along the patch boundaries
numElementsA = length(elementsA);
numElementsB = length(elementsB);

vertexA = zeros(1, length(elementsA));
vertexB = zeros(1, length(elementsB));

if (edgeA == 1) || (edgeA == 3) %horizontal edge
    vertIndexA = 3;
else
    vertIndexA = 4;
end

if (edgeB == 1) || (edgeB == 3) %horizontal edge
    vertIndexB = 3;
else
    vertIndexB = 4;
end

%loop over the elements in A
for elementIndex = 1:numElementsA    
    vertexA(elementIndex) = PHTelemA(elementsA(elementIndex)).vertex(vertIndexA);
end

%loop over the elements in B
for elementIndex = 1:numElementsB
    vertexB(elementIndex) = PHTelemB(elementsB(elementIndex)).vertex(vertIndexB);
end

newVertexB = setdiff(vertexA,vertexB);
newVertexA = setdiff(vertexB,vertexA);

vertexA = [0,vertexA];
vertexB = [0,vertexB];

%find the elements in A that need to be refined
for vertIndex = 1:length(newVertexA)
    for elementIndex = 1:numElementsA
        if (vertexA(elementIndex+1) > newVertexA(vertIndex)) && (vertexA(elementIndex) < newVertexA(vertIndex))
            elemRefA(elementsA(elementIndex)) = 1; %mark the corresponding element for refinement
            break
        end
    end
end

%find the elements and quads in B that need to be refined
for vertIndex = 1:length(newVertexB)
    for elementIndex = 1:numElementsB
        if (vertexB(elementIndex+1) > newVertexB(vertIndex)) && (vertexB(elementIndex)<newVertexB(vertIndex))            
            elemRefB(elementsB(elementIndex)) = 1; %mark the corresponding element for refinement
            break
        end
    end
end           
