function [image_out cp_imgIn] = resize_resample(export, im_dir, out_dir,out_fmt, image_in, scalefactor, intMethod)
    % This function is a wrapper that Ethan Duwell wrote for calling 
    % imresize to resize/resample images.
    % Input Variables:
    % -----------------
    % export: save an output image file in the output directory? (1 or 0)
    
    % im_dir: full path to directory where your input image file lives..
    
    % out_dir: full path to directory where you want to save the output image (if
    % desired..) if = "same" will set equal to im_dir
    %
    % out_fmt: desired output format/extension. If = "same", will assign to
    % be the same as input format/ext .. 
    %
    % image_in: filename of the input image+extension
    %
    % scalefactor: the scaling factor you want to apply (if <1 --> shrinks
    % image by factor, if >1 --> output image is bigger.. dimensions multiplied by factor..)
    %
    % intMethod: the interpolation method to use. ('nearest' for nearest neighbor is probably
    % preferable..)
    %
    %
    % Output Variables:
    % -----------------
    % image_out: output image matrix
    % NOTE: if export = 1 there will also be copies of the output image saved
    % in the specified output directory.
    
    %% Fetch The Input Image

    % save the starting directory location (so we can find our way back..)
    start_dir = pwd;
    
    % Go to the image directory..
    cd(im_dir)

    % Read in the input image..
    image_in_str = image_in; % save the name as a string first.. 
    image_in = imread(image_in);
    cp_imgIn = image_in; % assign to the output variable...
    
    %% Apply the requested rescaling/resizing using imresize()
    image_out = imresize(image_in,scalefactor,intMethod);

    %% If requested, save output image copy in output directory
    if export == 1
        % go to the output directory ..

        % if out_dir == "same", reassign to equal im_dir
        if out_dir == "same"
            out_dir = im_dir;
        end
        cd(out_dir);

        % create descriptive output image name specifying scaling factor
        % and the interpolation method used .. 

        % a. extract input filename root string from the extension 
            [filepath,fname,ext] = fileparts(image_in_str); % seperate parts using fileparts

            % cut the "." off the extension
            ext = extractAfter(ext,".");
        % b. get scalefactor as a string. Replace "."s in number with "_" if
        %    present..
            scalefactor_str = num2str(scalefactor);
            scalefactor_str = strrep(scalefactor_str,'.','_');
            
        % if out_fmt == "same", reassign to equal input image ext
            if out_fmt == "same"
                out_fmt = ext;
            end

        % c. contatenate parts extracted above along with the interp method 
        % and format extension to form the descriptive output filename..
            out_name = strcat(fname,"_rs",scalefactor_str,"_",intMethod,".",out_fmt);

        % d. write the file..        
            imwrite(image_out,out_name,out_fmt);
    end

    % go back to where you came from ..
    cd(start_dir);
end

