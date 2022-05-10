function globdat = solveStep(props,globdat)

switch props.nlSolver.type
    
    case 'staggered'
        
        globdat = staggeredSolver(props,globdat);

    case 'acclStaggered'
        
        globdat = acclStaggeredSolver(props,globdat);
        
    case 'nonlin'
        
        globdat = nonlinSolver(props,globdat);
        
    case 'qNonlin'
        
        error('Not yet implemented!')

    case 'nonlinArcLength'

        globdat = nonlinArcLength(props,globdat);
        
    otherwise
        
        error('Wrong option!')
        
end

end

