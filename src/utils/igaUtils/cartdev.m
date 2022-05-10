function[Basis]=cartdev(PHTelem,controlPts,geometry)
% Compute the derivatives, shape function, volume of each element and the gauss cordinates for each element.
% Uses PHTelem structure

p = geometry.p;
q = geometry.q;
ngaussX = geometry.ngaussX;
ngaussY = geometry.ngaussY;
numPatches = geometry.numPatches;

[GptsX,GWtsX]=GaussQuad(ngaussX);
[GptsY,GWtsY]=GaussQuad(ngaussY);

% 1D bernstein polynomials evaluated at the Gauss points on the master element
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

% Initialize the volume, dgdx and shape arrays
Basis = struct;
Basis.volume = cell(1,numPatches);
Basis.dgdx = cell(1,numPatches);
Basis.shape = cell(1,numPatches);
Basis.gaussCord = cell(1,numPatches);

for indexPatch = 1:length(PHTelem)
    
    Basis.volume(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.dgdx(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.shape(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    Basis.gaussCord(indexPatch) = {cell(length(PHTelem{indexPatch}),1)};
    
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            xmin = PHTelem{indexPatch}(i).vertex(1);
            xmax = PHTelem{indexPatch}(i).vertex(3);
            ymin = PHTelem{indexPatch}(i).vertex(2);
            ymax = PHTelem{indexPatch}(i).vertex(4);
            
            % The jacobian of the transformation from [-1,1]x[-1,1] to [xmin, xmax]x [ymin, ymax]
            scalefac = (xmax - xmin)*(ymax - ymin)/4;
            nument = size(PHTelem{indexPatch}(i).C,1);
            nodes = PHTelem{indexPatch}(i).nodes(1:nument);
            cpts = controlPts{indexPatch}(nodes,1:2);
            wgts = controlPts{indexPatch}(nodes,3);
            kgauss = 0;
            for igaussX=1:ngaussX
                for igaussY=1:ngaussY
                    kgauss = kgauss + 1;
                    dRdx = (PHTelem{indexPatch}(i).C)*squeeze(dBdu(igaussX,igaussY,:))*2/(xmax-xmin);
                    dRdy = (PHTelem{indexPatch}(i).C)*squeeze(dBdv(igaussX,igaussY,:))*2/(ymax-ymin);
                    RR = (PHTelem{indexPatch}(i).C)*squeeze(R(igaussX,igaussY,:));
                    RR = RR .* wgts;
                    dRdx = dRdx .* wgts;
                    dRdy = dRdy .* wgts;
                    
                    w_sum = sum(RR);
                    dw_xi = sum(dRdx);
                    dw_eta = sum(dRdy);
                    
                    dRdx = dRdx/w_sum - RR*dw_xi/w_sum^2;
                    dRdy = dRdy/w_sum - RR*dw_eta/w_sum^2;
                    RR = RR/w_sum;
                    
                    % multiply by the jacobian of the transformation from reference
                    % space to the parameter space
                    dR  = [dRdx';dRdy'];
                    dxdxi = dR*cpts;
                    coord = RR'*cpts;
                    Basis.gaussCord{indexPatch}{i}(kgauss,1:2) = coord(:);
                    % Solve for first derivatives in global coordinates
                    dR = dxdxi\dR;
                    J = abs(det(dxdxi));
                    
                    Basis.volume{indexPatch}{i}(kgauss) = J*scalefac*GWtsX(igaussX).*GWtsY(igaussY);
                    Basis.shape{indexPatch}{i}(kgauss,:) = RR(:);
                    Basis.dgdx{indexPatch}{i}(kgauss,1:2,:) = dR(1:2,:);
                    Basis.gaussCord{indexPatch}{i}(kgauss,3) = J*scalefac*GWtsX(igaussX).*GWtsY(igaussY);
                end % end igaussY
            end % end igaussX
        end
    end
end

