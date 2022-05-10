function globdat = checkExit(props,globdat)

if globdat.ts.step == props.tStep.nsteps || ...
   globdat.lodi(globdat.ts.step+1,1)/max(globdat.lodi(:,1)) < 0.02

% Set globdat.active to false to exit time-stepping
globdat.active = false;

end

end

