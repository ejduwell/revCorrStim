function RevCorFixImGen(window, backwindow_fix,backwindow,backwindow_cap, fixrect,occlines,penW,grayblack,crop,cropsz,method,gray,mainsrc, maindest,baseOutDir,out_fmt)
%Generate image of fixation only condition

% Clear the windows..
Screen('FillRect', window, gray);  % blank out the onscreen window
%Screen('FillRect', backwindow, gray);  % blank out the offscreen window
Screen('FillRect', backwindow_fix, gray);  % blank out the offscreen window
Screen('FillRect', backwindow_cap, gray);  % blank out the offscreen window

%Draw Fixation Markers
[backwindow_fix,backwindow] = PORCfix4a(backwindow_fix,backwindow,fixrect,occlines,penW,grayblack);

% Write to main window/rotate..
Screen('DrawTexture', backwindow_cap, backwindow, mainsrc, maindest, 45, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
% Screen('Flip', window);% show linking display %

% Capture window as image..
img_snap =Screen('GetImage', backwindow_cap); %capture image of backwindow as it stands now...
% Convert to greyscale
img_snap = rgb2gray(img_snap);


% Save image to file
% contatenate parts extracted above along with the interp method 
% and format extension to form the descriptive output filename..
out_name = strcat(baseOutDir,"/fixation.png");

if crop == 1
    img_snap = imgSizeEdit(img_snap,cropsz,method);
end

imwrite(img_snap,out_name,out_fmt);

% Clear the windows..
%Screen('FillRect', window, gray);  % blank out the onscreen window
Screen('FillRect', backwindow, gray);  % blank out the offscreen window
Screen('FillRect', backwindow_cap, gray);  % blank out the offscreen window

end