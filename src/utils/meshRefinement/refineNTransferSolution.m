function globdat = refineNTransferSolution(props,globdat)

% Transfer solution, global to local
[globdat.state,globdat.state0,globdat.state00,globdat.Dstate] = ...
    TransferNodalDataGlob2Loc(globdat.state,globdat.state0, ...
    globdat.state00,globdat.Dstate,globdat.PHTelem,globdat.mesh,props.geom.dim);               

% Refine mesh and transfer solution, local to global
[globdat.PHTelem,globdat.ctrlPts,globdat.mesh,globdat.markRef,globdat.state,globdat.state0,globdat.state00, ...
    globdat.Dstate] = RefineTransferNodalDataGlob2Loc(globdat.markRef, ...
    props.geom,globdat.state,globdat.state0,globdat.state00,globdat.Dstate,globdat.ctrlPts,...
    globdat.PHTelem,globdat.mesh);
            
% Compute basis functions for the new mesh
globdat.Basis = cartdevRefinement(globdat.PHTelem,globdat.ctrlPts,props.geom,...
                                      globdat.Basis,globdat.markRef);

% Set new ndofs
globdat.ndofs = length(globdat.state);

end