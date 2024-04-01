function [MatchArray, NonMatchArray, MatchIndexs, NonMatchIndexs] = compareCellStrLists(cellArray1,cellArray2)
% Accepts two N by 1 cell arrays of strings: cellArray1 and cellArray2
% Each string in cellArray1 is checked against cellArray2 to see if it is
% present. If present, it is saved in MatchArray. If not, it is saved in
% NonMatchArray. Row indices of matches are saved in MatchIndexs. Row
% indices of non-matches are saved in NonMatchIndexs.

% Preallocate
MatchArray={};
NonMatchArray={};
MatchIndexs={};
NonMatchIndexs={};
for ii=1:size(cellArray1,1)
    fileString=cellArray1{ii,1};
    if ~any(strcmp(cellArray2,fileString))
        NonMatchArray{end+1,1}=fileString;
        NonMatchIndexs{end+1,1}=ii;
    else
        MatchArray{end+1,1}=fileString;
        MatchIndexs{end+1,1}=ii;
    end
end

end