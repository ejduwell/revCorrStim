% Clear the workspace and the screen
sca;
close all;
clear;

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create onscreen window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
screenrect = Screen('Rect',0);	% get the size of the display screen

%scale_f = 1;
scale_f = 0.5;

% ejd attempt to rescale whole double screen into a single screen..
%screenrect = screenrect * scale_f; %ORIG
screenrect(3) = screenrect(3) * scale_f;

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');
[window, windowRect] = PsychImaging('OpenWindow', 0, GrayIndex(0,0.5), screenrect);	% Open generic on-screen window

Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    % allow for transparency!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Make a base Rect of 200 by 250 pixels
baseRect = [0 0 200 500];

% For Ovals we set a miximum diameter up to which it is perfect for
maxDiameter = max(baseRect) * 1.01;

% Center the rectangle on the centre of the screen
centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);

% Set the color of the rect to red
rectColor = [1 0 0];

% Draw the rect to the screen
%Screen('FrameOval', window, [], centeredRect, maxDiameter);
Screen('FrameOval', window, 0, centeredRect,2);

% Flip to the screen
Screen('Flip', window);

% Wait for a key press
KbStrokeWait;

% Clear the screen
sca;