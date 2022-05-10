%-------------------------------------------------------------------------%
% Sol = ApplyDBC(Sol,time,constraints) applies time-dependent Dirichlet 
% Boundary Condition to the solution vector 'Sol'. 
%
% INPUT:  Sol            -> Solution vector
%         time           -> Current time
%         constraints    -> Matrix containing [DirichletDOFS; Factor]
%                           Factor = 0 : homogeneous Dirichlet condition
%                           Factor = 1 : ihhomogeneous Dirichlet condition
%
% METHOD: Sol(Dirichlet DOFS) = time * Factor;
%
% OUTPUT: Sol            -> Solution vector with Dirichlet values applied
%
% AUTHOR: Ritukesh Bharali (ritukesh.bharali@chalmers.se)
%         Materials and Computational Mechanics,
%         Department of Industrial and Material Science,
%         Chalmers University of Technology, Gothenburg, Sweden.
%         
% DATE:   06.01.2021
%-------------------------------------------------------------------------%

function Sol = ApplyDBC(Sol,time,constraints)

Dirichlet_dofs           = constraints(:,1);
Dirichlet_factor         = constraints(:,2);
Sol(Dirichlet_dofs,1)      = time * Dirichlet_factor;

end