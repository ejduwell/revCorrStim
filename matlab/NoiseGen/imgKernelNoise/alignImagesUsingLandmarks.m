function alignedImage = alignImagesUsingLandmarks(image, referenceLandmarks, detectedLandmarks)
    % Compute the transformation required for alignment
    tform = fitgeotrans(detectedLandmarks, referenceLandmarks, 'affine');
    % Apply the transformation
    alignedImage = imwarp(image, tform, 'OutputView', imref2d(size(image)));
end