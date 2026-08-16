function mov = convert_videoToMat(video_path)
    % Load video and create cell array
    vidObj = VideoReader(video_path);

    disp(['Frame rate: ', num2str(vidObj.FrameRate), ' fps']);

    % Preallocate cell array
    numFrames = floor(vidObj.Duration * vidObj.FrameRate);
    mov = cell(1, numFrames);

    % Process each frame
    for k = 1:numFrames
        frame = read(vidObj, k);           % Read frame
        resizedFrame = imresize(frame, [112 112]); % Resize to 112x112
        mov{k} = resizedFrame;             % Store in cell array
    end

    % Save cell array to .mat file with the same name as the video file
    [videoDir, videoName, ~] = fileparts(video_path);
    save(fullfile(videoDir, [videoName, '.mat']), 'mov');
end
