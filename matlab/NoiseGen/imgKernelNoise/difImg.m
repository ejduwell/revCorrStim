function [difImg1,difImg2] = difImg(imgIn1,imgIn2,outTag)

im1=double(imread(imgIn1));
im2=double(imread(imgIn2));

difImg1=im2-im1;
difImg2=im1-im2;

difImg1=uint8(rescale(difImg1,0,255));
difImg2=uint8(rescale(difImg2,0,255));

figure
imshow(difImg1);
figure
imshow(difImg2);

imwrite(difImg1,strcat(outTag,"_1_R-L",".png"));
imwrite(difImg2,strcat(outTag,"_2_L-R",".png"));

end