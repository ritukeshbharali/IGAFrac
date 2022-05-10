function [quadRef] = refineMeshCoordinates(PHTelem,controlPts,quadList,p,q,xmin,xmax,ymin,ymax)
%marks for refinement all the elements on the PHTelem mesh which are inside the
%rectangle (xmin,xmax) x (ymin,ymax)

numPts = 2; %number of evaluation points to use on each edge
uref = linspace(-1,1,numPts);
vref = linspace(-1,1,numPts);

% 1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u] = bernstein_basis(uref,p);
[B_v] = bernstein_basis(vref,q);

R = zeros(numPts,numPts,(p+1)*(q+1));

basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        R(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
    end
end

numQuads = size(quadList,1);

% Initialize the marked patch array to zeros (i.e. no refinement)
quadRef = zeros(1, numQuads);


for i=1:numQuads
    for j=1:4
        % Determine the edges of the element in the physical space
        coord_south = zeros(numPts,2);
        coord_east = zeros(numPts,2);
        coord_north = zeros(numPts,2);
        coord_west = zeros(numPts,2);
        
        curElem = quadList(i,j);
        nument = min(size(PHTelem(curElem).C,1),length(PHTelem(i).nodes));
        nodes = PHTelem(curElem).nodes(1:nument);
        cpts = controlPts(nodes,1:2);
        wgts = controlPts(nodes,3);
        
        for jj=1:numPts
            % Compute the points on the element south edge
            N = squeeze(R(jj,1,:));
            NS = (PHTelem(curElem).C(1:nument,:))*N;
            NS = NS.*wgts;
            w_sum = sum(NS);
            NS = NS/w_sum;
            
            %compute the points on the element east edge
            N = squeeze(R(numPts,jj,:));
            NE = (PHTelem(curElem).C(1:nument,:))*N;
            NE = NE.*wgts;
            w_sum = sum(NE);
            NE = NE/w_sum;
            
            %compute the points on the element north edge
            N = squeeze(R(numPts-jj+1,numPts,:));
            NN = (PHTelem(curElem).C(1:nument,:))*N;
            NN = NN.*wgts;
            w_sum = sum(NN);
            NN = NN/w_sum;
            
            %compute the points on the element west edge
            N = squeeze(R(1,numPts-jj+1,:));
            NW = (PHTelem(curElem).C(1:nument,:))*N;
            NW = NW.*wgts;
            w_sum = sum(NW);
            NW = NW/w_sum;
            
            coord_south(jj,:) = NS'*cpts;
            coord_east(jj,:) = NE'*cpts;
            coord_north(jj,:) = NN'*cpts;
            coord_west(jj,:) = NW'*cpts;
        end
        coords = [coord_south; coord_east; coord_north; coord_west];
        
        %check if the coords are inside the marked region
        for iCoord = 1:size(coords,1)
            xCoord = coords(iCoord,1);
            yCoord = coords(iCoord,2);
            if (xmin <= xCoord) && (xCoord <= xmax) && (ymin <= yCoord) && (yCoord <= ymax)
                quadRef(i) = 1;
                break;
            end            
            
        end
        
    end
    
end

