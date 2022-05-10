% Function to have more scope to play around with the linear solver

function x = LinearSolver(A,b,LinSolver)

switch LinSolver.type
    
    case 'direct'
        
        x = A\b;
        
    case 'gmres'
        
        error('Only direct solver implemented!')
        
    otherwise
        
        error('Only direct solver implemented!')
        
end

end

