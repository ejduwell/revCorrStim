% RevCorr_main7
% -------------------------------------------------------------------------
% This is the main script for running reverse correlation object grouping
% experiments
%
% This script is intended to be run 'interactively.' To run it, simply call:
% 'RevCorr_main7' on the MATLAB command window. You will then be prompted
% to enter various parameters to run the experiment in the manner of your
% chosing..
%
% Written by: E.J. Duwell, PhD while working as a postdoc in 
% SNAP Lab run by Adam S. Greenberg's lab at MCW, Milwaukee WI
% 4/29/24
% -------------------------------------------------------------------------

clc; % clear command line interface..
%% Parameters

% get location of the "main" directory (directory housing this script..)
current_version = 'RevCorr_main7.m'; % will need to update with versioning..
mlab_dir = fileparts(which(current_version));
% cd to it..
cd(mlab_dir);

% get main dir too..
cd ..
cd ..
mainDir=pwd;
cd(mlab_dir);%return to mlab_dir


diaryfile=strcat(mlab_dir,"/commandWindowDiary");
if isfile(diaryfile)
    % if one exists, delete it..
    disp("detected old diary file..")
    disp("deleting old diary file and making a new one..")
    disp("starting the diary file..");
     delete(diaryfile);
     diary(diaryfile);
else
    disp("starting the diary file..");
     diary(diaryfile);
end

disp(" ");
disp(datetime);
disp(" ");
rng("shuffle")

%% Set up welcome banner

% Vector-o-rediculous-greetings
greeting_vect = ["Howdy, partner!","Ahoy, matey!","Avast m'harties!","Peek-a-boo!","Ello, gov'nor!","What's crackin?","Sup, homeslice?","Greetings and salutations.","Howdy, howdy, howdy!","Whaddup?!?","Shiver me timbers!","Ello, mate.","Why, hello there!","Aloha!", "Shalom","Que pasa Mufasa?","Bonjour!","Hallo", "Hei Hei! Hvorden har du det?","Ciao","Konnichiwa","Hiya!","Howdy-doody!","Yeeeeeeee Haaaawww!", "Yoooooouuu Hooooooo!!!!","Top of the mornin, to ya!","What's the word, hummingbird?","Hola, como estas?", "Hello! Is this the person to whom I'm speaking?", "What it do buckaroo?", "Greetings Comrade!","Good day mortals!","What's the word baby birds?","Guten Tag!","What's kickin, chicken?", "What's shakin, bacon?", "So... we meet at last!", "Greetings friend!","Peace be with you","Ayyyyyyoooooooooooooo!", "Shhaaaaaazzzzaaaaaammm!", "I BELIEVE IN A THING CALLED LOVE! JUSTLISTENTOTHERHYTHMOFYAHEEAAARRRT!!... ehhem.. excuse me...", "And I Said .. Heeyyyyaaaayahhyayaya..", "My regards to your superiors.", "Welcome to my humble abode.", "Greetings Earthlings!", "Fancy meeting you here!", "Welcome to chaos!", "Shalomie my homie!", "Wake up Neo... The Matrix has you...", "It was best of times. It was the worst of times.", "I'm glad that you are here with me. Here at the end of all things.", "Lets get down to business! To complete.. some runs!", "Haaaallp! Get me out of here!","Pax Hominibus.. Peace on the Minibus.","Stop monkeying around and get to work!", "This experiment bears the indelible stamp of our lowly origins..", "What is up my scallywagz?","Och aye the noo!","Velkommen!","Hou ar ye?", "I came in like a wreeeeecking balll!!!"];
n_greetings = length(greeting_vect);
ran_idx = randi([1,n_greetings]);

greeting = greeting_vect(ran_idx);
brk_line = '-------------------------------------------------------------------------------------------------------------------\n';
greeting2 = '*********************** This is the main script for running reverse correlation experiments ***********************';
ref_str = greeting2;
g_len = strlength(greeting);
g_len2 = strlength(ref_str);
d_len = g_len2-g_len;
if ~mod(d_len,2) == 1
    len_add1 = round(d_len/2)-1;
    len_add2 = len_add1;
else
    len_add1 = round(d_len/2)-2;
    len_add2 = round(d_len/2)-1;
end
chk = "-";
chk_init1 = chk;
for ii = 1:(len_add1-1)
chk_init1 = strcat(chk_init1,chk);
end
clear ii
chk_init2 = chk;
for ii = 1:(len_add2-1)
chk_init2 = strcat(chk_init2,chk);
end
clear ii

greeting = strcat(chk_init1," ",greeting," ",chk_init2,'\n');
greeting1 = greeting;


%% Run interactive CLI for selecting/starting the experiment 

% 1) Welcome Screen
pauseDuration=0.05; % Controls the rate at which the banner image "scrolls out"
asciiMatPath=strcat(mlab_dir,"/asciiBannerIms/asciiArt.mat"); % path to mat of asciis..
fprintf('\n\n')
randAsciiArt(asciiMatPath,pauseDuration)
fprintf('\n')
fprintf(greeting1)
fprintf(greeting2)
fprintf('\n')
versionmsg = strcat("(",current_version, ", Written by: E.J. Duwell, PhD)");
versionmsg = cln_cmdlinemsg(ref_str,versionmsg," ");

instr_pgMsg = "Read: https://github.com/ejduwell/revCorrStim/blob/main/README.md";

instr_pgMsg = cln_cmdlinemsg(ref_str,instr_pgMsg," ");
fprintf(versionmsg)
fprintf(instr_pgMsg)
fprintf(brk_line)
continuemsg = cln_cmdlinemsg(ref_str,"Press Enter/Return to Continue..."," ");
fprintf(continuemsg)
input("",'s')

% 2) Pick a Version
fprintf('\n\n')

%RevCorr_QSTmain6
RevCorr_QSTmain7

diary off;
clear;

%% Functions
function  [output] = cln_cmdlinemsg(longest_text,input,filler)
    greeting2 = longest_text ;
    greeting = input;
    g_len = strlength(greeting);
    g_len2 = strlength(greeting2);
    d_len = g_len2-g_len;
    if ~mod(d_len,2) == 1
        len_add1 = round(d_len/2);
        len_add2 = len_add1;
    else
        len_add1 = round(d_len/2)-1;
        len_add2 = round(d_len/2);
    end
    chk = filler;
    chk_init1 = chk;
    for ii = 1:(len_add1)
    chk_init1 = strcat(chk_init1,chk);
    end
    clear ii
    chk_init2 = chk;
    for ii = 1:(len_add2)
    chk_init2 = strcat(chk_init2,chk);
    end
    clear ii

    greeting = strcat(chk_init1," ",greeting," ",chk_init2,'\n');
    output = greeting;
end

function  [output] = add_space(longest_text,input,filler)
    greeting2 = longest_text ;
    greeting = input;
    g_len = strlength(greeting);
    g_len2 = strlength(greeting2);
    d_len = g_len2-g_len;

    len_add1 = round((g_len2/4)-6);
        
    chk = filler;
    chk_init1 = chk;
    for ii = 1:(len_add1-1)
    chk_init1 = strcat(chk_init1,chk);
    end
    clear ii

    greeting = strcat(chk_init1,greeting," ");
    output = greeting;
end
