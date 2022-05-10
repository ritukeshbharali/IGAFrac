function [dirichlet] = DirichletConstraints(PHTelem,geometry,NN)

% Some useful variables
p = geometry.p;
q = geometry.q;

% Define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

% TO-DO: Pre-allocate and then remove empty entries for speed up as the
% array re-sizing happens only once and not for every element on the
% boundary.

bottomEdge = [];
topEdge = [];
leftEdge = [];
rightEdge = [];

% Set the bottom boundary degrees of freedom
for patchIndex = 1:2
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_down)
                bottomEdge = [bottomEdge, PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];
            end
        end
    end
end

% Set the top boundary degrees of freedom
for patchIndex = 3:4
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_up)
                topEdge = [topEdge, PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
            end
        end
    end
end

% Set the left boundary degrees of freedom
for patchIndex = [1,4]
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_left)
                leftEdge = [leftEdge, PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
            end
        end
    end
end

% Set the right boundary degrees of freedom
for patchIndex = [2,3]
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_right)
                rightEdge = [rightEdge, PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
            end
        end
    end
end

% Eliminate duplicate entries
Left = unique(leftEdge,'stable');
Right = unique(rightEdge,'stable');
Top = unique(topEdge,'stable');
Bottom = unique(bottomEdge,'stable');

% SOMDATTA: Dirichlet dofs
% dofs = [2*Bottom, ...
%         2*Top, ...
%         2*Left-1];
%               
% % Multiplying factor for dirichlet dofs (0 - homogeneous, 1 - otherwise)              
% vals = zeros(1,length(dofs));
% l_bottom = length(Bottom);
% l_top = length(Top);
% vals(l_bottom+1:l_bottom+l_top) = 1.0;

% RITUKESH: I changed the constraints to set I used in COMSOL
% dofs = [2*Bottom-1, 2*Bottom ...
%         2*Top-1, 2*Top ...
%         2*Left, ...
%         2*Right];

% dofs = [2*Bottom-1, 2*Bottom ...
%         2*Top-1, 2*Top ...
%         2*Left ];    

dofs = [2*Bottom-1, 2*Bottom ...
        2*Top ];   

vals     = zeros(1,length(dofs));
l_bottom = length(Bottom);
l_top    = length(Top);
vals(2*l_bottom+1:2*l_bottom+l_top) = 1.0;

% Fill in the dirichlet dofs and their values
dirichlet = [dofs; vals]';

end