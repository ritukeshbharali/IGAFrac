function[shape,dgdx,volume,gaussCord] = CalcDerv(i,PHTelem,controlPts,ngaussX,ngaussY,GWtsX,GWtsY,R,dBdu,dBdv)             

xmin = PHTelem(i).vertex(1);
xmax = PHTelem(i).vertex(3);
ymin = PHTelem(i).vertex(2);
ymax = PHTelem(i).vertex(4);

% The jacobian of the transformation from [-1,1]x[-1,1] to [xmin, xmax]x [ymin, ymax]
scalefac = (xmax - xmin)*(ymax - ymin)/4;
nument = size(PHTelem(i).C,1);
nodes = PHTelem(i).nodes(1:nument);
cpts = controlPts(nodes,1:2);
wgts = controlPts(nodes,3);
kgauss = 0;
shape = zeros(ngaussX*ngaussY,nument);
gaussCord = zeros(ngaussX*ngaussY,2);
volume = zeros(ngaussX*ngaussY);
dgdx = zeros(ngaussX*ngaussY,2,nument);
for igaussX=1:ngaussX
    for igaussY=1:ngaussY
        kgauss = kgauss + 1;
        dRdx = (PHTelem(i).C)*squeeze(dBdu(igaussX,igaussY,:))*2/(xmax-xmin);
        dRdy = (PHTelem(i).C)*squeeze(dBdv(igaussX,igaussY,:))*2/(ymax-ymin);
        RR = (PHTelem(i).C)*squeeze(R(igaussX,igaussY,:));
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
        gaussCord(kgauss,:) = coord(:);
        % Solve for first derivatives in global coordinates
        dR = dxdxi\dR;
        J = abs(det(dxdxi));
        
        volume(kgauss) = J*scalefac*GWtsX(igaussX).*GWtsY(igaussY);
        shape(kgauss,:) = RR(:);
        dgdx(kgauss,1:2,:) = dR(1:2,:);
    end % end igaussY
end % end igaussX