function squareCrop(fullImPath,scalefactor)

[filePath, fileName, fileExt] = fileparts(fullImPath);

image_in=strcat(fileName,fileExt);
[image_out, cp_imgIn] = resize_resample(0, filePath, filePath,"same", image_in, scalefactor, 'nearest');

image_out = rgb2gray(image_out);
[xdim,ydim] = size(image_out);
desiredSize=min([xdim,ydim]);

[imgOut2] =  imgSizeEdit(image_out,desiredSize,'nearest');





out_name=strcat(filePath,"/",fileName,"_rs",num2str(scalefactor),"_sqr",fileExt);

imwrite(imgOut2,out_name,fileExt(2:end));

end