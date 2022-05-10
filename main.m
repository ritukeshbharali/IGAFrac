%=========================================================================%
% IGAFrac is a MATLAB software package for PHT spline-based adaptive 
%             Iso-Geometric Analysis (IGA) of fractures.
% 
% Release:     10.05.2022
% Authors:     Ritukesh Bharali (Chalmers University of Technology, Sweden)
%              Somdatta Goswami (Brown University, USA)
%=========================================================================%

function globdat = main(casename)

% If no casename is provided, ask for one
if nargin ~= 1
    casename = input('Enter problem to simulate: ','s');
end

% Clear workspace
close all; clc
format long;
warning('off')
tic

%-------------------------------------------------------------------------%
% ADD PATH TO DIRECTORIES
%-------------------------------------------------------------------------%

disp(' - Adding path to directories')

% Path to external libraries and source folders
addpath(genpath('./ext'));
addpath(genpath('./src/io'));
addpath(genpath('./src/materials'));
addpath(genpath('./src/solver'));
addpath(genpath('./src/utils'));

% Remove path to femodel and inputData to avoid function clashes
rmpath(genpath('./src/igamodel'))
rmpath(genpath('./inputData'))

% Path to input files
addpath(fullfile('./inputData/',casename));

% Create an output directory
dir_output = sprintf('./output/%s', casename);
if ~exist(dir_output, 'dir')
    mkdir(dir_output)
end

% Remove path to output, and add only case specific directory
rmpath(genpath(dir_output))
addpath(fullfile(dir_output,casename))

%-------------------------------------------------------------------------%
% SETUP RUNTIME PARAMETERS, MESH
%-------------------------------------------------------------------------%

disp(' - Reading user input and setting up runtime parameters and mesh')
props          = inputData;

% Add path to the right femodel
addpath(fullfile('src/igamodel/',props.femodel.type))

% Initialise the global database
globdat        = setupProblem(props);

%-------------------------------------------------------------------------%
% TIME-STEPPING
%-------------------------------------------------------------------------%

while globdat.active
           
    % Solve a step
    globdat = solveStep(props,globdat);
    
    % Check for exit
    globdat = checkExit(props,globdat);
    
    % Post-process
    globdat = postProcessStep(dir_output,props,globdat);
    
end

%-------------------------------------------------------------------------%
% SHUTDOWN
%-------------------------------------------------------------------------%

shutdownProblem(dir_output,props,globdat);
disp(' ')
disp([' - *** TIME ELAPSED = ',num2str(toc),' seconds.'])

%=========================================================================%
% End of main function
%=========================================================================%
end
