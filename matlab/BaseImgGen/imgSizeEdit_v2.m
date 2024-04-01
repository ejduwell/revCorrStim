function [imgOut] =  imgSizeEdit_v2(imgIn,desiredSize,method)

% --- Crop --- %
imgIn = imcrop(imgIn,[(size(imgIn,2)-size(imgIn,1))/2, 0, size(imgIn,1), size(imgIn,1)]);

% --- Make sure dims are the same --- %
if size(imgIn ,1) ~= size(imgIn ,2)
    imgSize   = size(imgIn);
    [val,idx] = max(imgSize);
    switch idx
        case 1
            imgIn (val,:) = [];
        case 2
            imgIn(:,val) = [];
    end
end
% --- Resize & Store --- %
%imgOut = imresize(imgIn ,[desiredSize, desiredSize],method);

imgOut = imresize(imgIn ,[desiredSize, desiredSize],method,Antialiasing=false);

end