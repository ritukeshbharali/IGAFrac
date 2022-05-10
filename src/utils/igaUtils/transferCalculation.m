function[shape,dgdx,Refelem] = transferCalculation(ElmtList,PHTelem,controlPts,p,q,xmin,xmax,ymin,ymax,nGauss)

dim = 2;
if nargin<10
    nGauss = p+1;
end
ngaussX = nGauss;
ngaussY = nGauss;
[Gpts,~]=GaussQuad(ngaussX);
toleq = 1e-10;
counter = 1;
for iCount = 1:length(ElmtList)
    if isempty(PHTelem(ElmtList(iCount)).children)% This loop is for the elements that have been refined in the other mesh and has the reference elemnt number affected in the current mesh
        leaf(counter) = ElmtList(iCount);% Leaf contains the elemnt number ogf the most refined element. The elemt number stored in leaf does not have children
        vertex(counter,:) = PHTelem(ElmtList(iCount)).vertex;
        counter = counter + 1;
    else
        leafTemp = PHTelem(ElmtList(iCount)).children;% Children have to be checked since in the current miter, there might be refinement in other Mesh
        for  ileaf = 1 : 4 % Coz each element has 4 children
            leaf(counter) = leafTemp(ileaf);
            vertex(counter,:) = PHTelem(leafTemp(ileaf)).vertex;
            counter = counter + 1;
        end
    end
end
kgauss=0;
shape = zeros(ngaussX*ngaussY,(p+1)*(q+1));
dgdx = zeros(ngaussX*ngaussY,dim,(p+1)*(q+1));
Refelem = zeros(1,ngaussX*ngaussY);
for ii=1:ngaussX
    for jj=1:ngaussY
        kgauss=kgauss+1;
        u_hat = Gpts(ii);
        v_hat = Gpts(jj);
        
        % Mapping (u_hat,v_hat) on [-1,1]x[-1,1] to (xi,eta) on[xmin,xmax]x[ymin,ymax]
        xi = u_hat*(xmax-xmin)/2+(xmax+xmin)/2;
        eta = v_hat*(ymax-ymin)/2+(ymax+ymin)/2;
        
        for iElmt = 1:counter-1
            if (xi <= vertex(iElmt,3)+toleq) && (xi >= vertex(iElmt,1)-toleq) && (eta <= vertex(iElmt,4)+toleq) && (eta >= vertex(iElmt,2)-toleq)
                elmtNum = leaf(iElmt);
                Pxmin = vertex(iElmt,1);
                Pxmax = vertex(iElmt,3);
                Pymin = vertex(iElmt,2);
                Pymax = vertex(iElmt,4);
                break;
            end
        end
        uu_hat = (2*xi - Pxmin - Pxmax)/(Pxmax-Pxmin);
        vv_hat = (2*eta - Pymin - Pymax)/(Pymax-Pymin);
        
        nument = size(PHTelem(elmtNum).C,1);
        nodes = PHTelem(elmtNum).nodes(1:nument);
        cpts = controlPts(nodes,1:2);
        wgts = controlPts(nodes,3);
        shape(kgauss,:) = phtBasis(uu_hat,vv_hat,PHTelem(elmtNum).C,p,q);
        [~,dR] = phtBasisIso(uu_hat,vv_hat,PHTelem(elmtNum).C,p,q,wgts);
        dR(1,:) = dR(1,:)*2/(Pxmax-Pxmin);
        dR(2,:) = dR(2,:)*2/(Pymax-Pymin);
        dxdxi = dR*cpts;
        dR = dxdxi\dR;
        dgdx(kgauss,:,:) = dR;
        Refelem(kgauss) = elmtNum;
    end %igaus
end %jgaus