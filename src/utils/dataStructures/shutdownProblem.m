function shutdownProblem(dir_output,props,globdat)

% Re-size lodi array and export to txt file
globdat.lodi(globdat.ts.step+2:end,:) = [];
writematrix(globdat.lodi,fullfile(dir_output,'lodi.txt'))  

% Export error vectors to txt file
for i = 1:length(globdat.errVec)
    writematrix(globdat.errVec{i},fullfile(dir_output,sprintf('errVec%g', i)))  
end

% Zip up the output folder
% zip([dir_output,'.zip'],dir_output)

end