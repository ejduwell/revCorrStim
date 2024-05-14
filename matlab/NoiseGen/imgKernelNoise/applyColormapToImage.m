function colorImg = applyColormapToImage(grayImg, colormapName)

    jnk=figure;
    % Ensure the grayscale image is in double format and normalized
    if ~isa(grayImg, 'double')
        grayImg = double(grayImg);
    end
    grayImg = (grayImg - min(grayImg(:))) / (max(grayImg(:)) - min(grayImg(:)));

    % Apply the specified colormap
    cmap = colormap(colormapName); % Get the specified colormap

    % Convert the normalized grayscale image to an index image
    %[~,~,~,idx] = gray2ind(grayImg, size(cmap, 1));
    %[~,~,~,idx] = gray2ind(grayImg, size(cmap, 1));
    [idx,cmap2] = gray2ind(grayImg,size(cmap, 1));

    % Convert the indexed image to a true color image using the colormap
    colorImg = ind2rgb(idx, cmap);

    % Ensure output is in double format
    colorImg = im2double(colorImg);

    %close any figure windows which might be created above..
    close(jnk);
end