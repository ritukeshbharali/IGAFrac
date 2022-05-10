function makeMovie(dir_output,filename,frameRate,step)

%filePattern = fullfile(locationFolder,'Step*.png');% Get a directory listing
%imageFiles = dir(filePattern);

writerObj = VideoWriter(fullfile(dir_output,[filename,'.avi']));% Open the video writer object.
writerObj.FrameRate = frameRate;
writerObj.Quality   = 95;
open(writerObj);

for iframe = 1 : step
    filename = sprintf('Step%d.png',iframe);
    frame = imread(filename);
    writeVideo(writerObj,frame);% Write this frame out to the AVI file
end

close(writerObj);% Close down the video writer object to finish the file

end