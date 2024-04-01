function landmarks = detectFacialLandmarks(img)
    % Initialize the face, eyes, nose, and mouth detectors
    faceDetector = vision.CascadeObjectDetector('FrontalFaceCART');
    eyeDetector = vision.CascadeObjectDetector('EyePairBig');
    noseDetector = vision.CascadeObjectDetector('Nose');
    mouthDetector = vision.CascadeObjectDetector('Mouth');

    % Detect the face
    bboxFace = step(faceDetector, img);
    if isempty(bboxFace)
        landmarks = [];
        return;
    end

    % Consider the first detected face
    faceRegion = imcrop(img, bboxFace(1, :));

    % Detect eyes, nose, and mouth within the face region
    bboxEyes = step(eyeDetector, faceRegion);
    bboxNose = step(noseDetector, faceRegion);
    bboxMouth = step(mouthDetector, faceRegion);

    % Calculate the center points of eyes, nose, and mouth
    landmarks = zeros(3, 2); % Initializing the landmarks array
    if ~isempty(bboxEyes)
        landmarks(1, :) = [bboxEyes(1) + bboxEyes(3)/2, bboxEyes(2) + bboxEyes(4)/2] + bboxFace(1, 1:2);
    end
    if ~isempty(bboxNose)
        landmarks(2, :) = [bboxNose(1) + bboxNose(3)/2, bboxNose(2) + bboxNose(4)/2] + bboxFace(1, 1:2);
    end
    if ~isempty(bboxMouth)
        landmarks(3, :) = [bboxMouth(1) + bboxMouth(3)/2, bboxMouth(2) + bboxMouth(4)/2] + bboxFace(1, 1:2);
    end
end
