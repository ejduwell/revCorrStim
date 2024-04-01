function struct_out = evalStrVect(strVect)
%evalStrVect Sequentially unpacks and executes the strings contained in 
% strVect as MATLAB code..
%   Iterates through the list of strings in the input "strVect" and
%   executes each as a MATLAB command using the "eval" function.

nreps = length(strVect);

for ii = 1:nreps
    cmd = strVect(ii);
    eval(cmd);
end

pause = "";

end