function [controlPts,PHTelem,dimBasis,quadList] = initElements(nurbs,p,q,numElemU,numElemV)

    p_init = nurbs.order(1)-1;
    q_init = nurbs.order(2)-1;
    % Refine into numberElemU by numberElemV knotspans
    knotU = linspace(0,1,numElemU+1);
    knotV = linspace(0,1,numElemV+1);
    
    numberElementsU = length(unique(knotU))-1;
    numberElementsV = length(unique(knotV))-1;
    
    % Increase polynomial order
    nurbs = nrbdegelev(nurbs,[p-p_init,q-q_init]);
    nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});
    % Repeat the knots to get C1 continuity
    nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});
    [controlPts,PHTelem,dimBasis,quadList] = genControlPtsNoRep(nurbs,p,q,numberElementsU,numberElementsV);
    
    [quadRef] = refineMeshCoordinates(PHTelem,controlPts,quadList,p,q,3.55+0.01,4.55-0.01,0.95+0.01,1.105-0.01);
    [quadList, PHTelem, controlPts,dimBasis] = refineMeshIso(quadRef,quadList,PHTelem,controlPts,p,q,dimBasis);

    [quadRef] = refineMeshCoordinates(PHTelem,controlPts,quadList,p,q,3.525+0.001,4.475-0.001,0.975+0.001,1.0525-0.001);
    [quadList,PHTelem, controlPts,dimBasis] = refineMeshIso(quadRef,quadList,PHTelem,controlPts,p,q,dimBasis);
      
end