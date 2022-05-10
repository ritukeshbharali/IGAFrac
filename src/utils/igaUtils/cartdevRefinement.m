function[Basis]=cartdevRefinement(PHTelem,controlPts,geometry,Basis,markRef)
% This evaluates the information only of the new elements formed

p = geometry.p;
q = geometry.q;
ngaussX = geometry.ngaussX;
ngaussY = geometry.ngaussY;
numPatches = geometry.numPatches;

[GptsX,GWtsX]=GaussQuad(ngaussX);
[GptsY,GWtsY]=GaussQuad(ngaussY);

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u,dB_u] = bernstein_basis(GptsX,p);
[B_v,dB_v] = bernstein_basis(GptsY,q);

dBdu = zeros(ngaussX, ngaussY, (p+1)*(q+1));
dBdv = zeros(ngaussX, ngaussY, (p+1)*(q+1));
R = zeros(ngaussX, ngaussY, (p+1)*(q+1));

% The derivatives of the 2D Bernstein polynomials at Gauss points on the master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
        R(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
    end
end

for indexPatch = 1:numPatches
    for indexRef = 1:length(markRef{indexPatch})
        if markRef{indexPatch}(indexRef)
            if isempty(PHTelem{indexPatch}(indexRef).children)
                [Basis.shape{indexPatch}{indexRef},Basis.dgdx{indexPatch}{indexRef},...
                    Basis.volume{indexPatch}{indexRef},Basis.gaussCord{indexPatch}{indexRef}] = ...
                    CalcDerv(indexRef,PHTelem{indexPatch},controlPts{indexPatch},ngaussX,ngaussY,GWtsX,GWtsY,R,dBdu,dBdv);               
            else
                childElmt = PHTelem{indexPatch}(indexRef).children;
                for ichild = 1:4
                    i = childElmt(ichild);
                    [Basis.shape{indexPatch}{i},Basis.dgdx{indexPatch}{i},Basis.volume{indexPatch}{i},...
                        Basis.gaussCord{indexPatch}{i}] = CalcDerv(i,PHTelem{indexPatch},controlPts{indexPatch},ngaussX,ngaussY,GWtsX,GWtsY,R,dBdu,dBdv);               
                end
                Basis.shape{indexPatch}{indexRef} = [];% Deactivating the parent element which has been refined
                Basis.dgdx{indexPatch}{indexRef} = [];
                Basis.volume{indexPatch}{indexRef} = [];
                Basis.gaussCord{indexPatch}{indexRef} = [];
            end
        end
    end
end
