
function makeCheckerboard()


% Clear the workspace and the screen
sca;
close all;
clear;

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create onscreen window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
screenrect = Screen('Rect',0);	% get the size of the display screen

scale_f = 0.5;

% ejd attempt to rescale whole double screen into a single screen..
%screenrect = screenrect * scale_f; %ORIG
screenrect(3)=1680;

%screenrect(3) = screenrect(3) * scale_f;

% Open an on screen window
screenNumber=0;
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, GrayIndex(0,0.5), screenrect);	% Open generic on-screen window

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white / 2;

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Define a simple 2*(chkrBrdSz) by 2*(chkrBrdSz) checker board, 
% If chkrBrdSz is set to 2 below, this will effectively be a four by
% four pixel image that we scale using onboard PTB filtering.
% if chkrBrdSz is set to 1, it will be a 2 by 2
% if chkrBrdSz is set to 3 it will be a 6 by 6... etc..

% chkrBrdSz controls the size of the checkerboard array 
% 
% (ie this x2 will be the width/height of the array im units of total # of checks)
chkrBrdSz=10; 
sizeArrayOut=10000; % square dimension of output array..

checkerboard = repmat(eye(2), chkrBrdSz, chkrBrdSz);

% Make the checkerboard into a texure (4 x 4 pixels)
checkerTexture = Screen('MakeTexture', window, checkerboard);

% We will scale our texure up to 90 times its current size be defining a
% larger screen destination rectangle\

% checkerRectScaler is a scaling factor that controls size of 
% checker array in terms of height/width of each checker in the array 
% (in pixels). 
checkerRectScaler=20; 


[s1, s2] = size(checkerboard);
dstRect = [0 0 s1 s2] .* checkerRectScaler;
dstRect = CenterRectOnPointd(dstRect, xCenter, yCenter);

% Draw the checkerboard texture to the screen. By default bilinear
% filtering is used. For this example we don't want that, we want nearest
% neighbour so we change the filter mode to zero
filterMode = 0;
Screen('DrawTextures', window, checkerTexture, [],dstRect, [], filterMode);

% Flip to the screen
Screen('Flip', window);

% Wait for a key press
KbStrokeWait;

% Clear the screen
sca;



end