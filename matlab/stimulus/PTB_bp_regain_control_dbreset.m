%function PTB_bp_regain_control_dbreset()
% Ethan Duwell wrote this function (March 2022) to call the necessary commands to regain
% control of mouse, keyboard and screen when debugging psych toolbox
% scripts/code/stimuli etc.. 
%
% Intended use: Ethan puts PTB_bp_regain_control after the "Catch"
% statement in a Try-Catch block spanning the entire script which is being
% debugged. That way, when something breaks the code, you will be bumped
% out into "Catch" and these commands will be called to give you back
% control of you computer. It also re-prints the last error message to help
% zero in on what is broken.
%
% Alternatively: This function can also be placed in the line preceeding a
% breakpoint while debugging. This will ensure that when the breakpoint is
% reached, you will be able to use your mouse, keyboard and see your Matlab
% command prompt/IDE window on screen..

db_dummy_file_path = '/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/data/';
db_dummy_file = 'XXX1.mat';
ListenChar(0);
ShowCursor;
Screen('CloseAll');
sca
path_delete(db_dummy_file_path,db_dummy_file); %delete dummy file during debug
clear all
close all
rethrow(lasterror); % Spit the last error message out on the command line so you can know what is broken
%end