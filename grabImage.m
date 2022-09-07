function ImBG_autoCont = grabImage(RCheckBox, GCheckBox, BCheckBox)
    [err, M, w, h] = uEye_camera(1); % fetch image
        if err>0
            error('Capturing an image failed!')
        end
    M(M<0)=M(M<0)+256; %Restructure incoming pixel values to fall between 0 and 256
    M = reshape(M, 3, w*h); %Create 3D matrix for 2D image in R, G and B
    I = reshape(M, 3, w, h);
    I = permute(I, [3 2 1]);
    RVal = RCheckBox.Value; % Get checkbox values for which colour channels to use
    GVal = GCheckBox.Value;
    BVal = BCheckBox.Value;
    Im_mono = 0.55*RVal*I(:,:,2)+0.35*GVal*I(:,:,3)+0.1*BVal*I(:,:,1); % Create monochrome image from the 3 colour channels
    ImBG_autoCont = rescale(Im_mono,0,256); % Rescale pixel values to maximise contrast, keep floating point precision
end