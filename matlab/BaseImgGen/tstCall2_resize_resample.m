% test pars ..
export = 1;
%im_dir = "/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir";
im_dir = "/Users/ggurariy/Library/CloudStorage/OneDrive-mcw.edu/ReverseCorrelation/BaseImgGen/texImgs";
out_dir = "same";
out_fmt = "same";
image_in = "gravel_highres.jpg";
scalefactor = 10;
intMethod = "nearest";

% test call .. 
[test_out, test_orignal] = resize_resample(export, im_dir, out_dir,out_fmt, image_in, scalefactor, intMethod);
figure
imshow (test_out)
figure
imshow (test_orignal)

pause = "";