function [itrParsCom_matches] = combvecLogicalIdx(itrPars,itParNamz,Eqlz,itrParsCom)


%make numeric indexes for the names in itParNamz
itParNamz_idx=linspace(1,length(itParNamz),length(itParNamz));

%preallocate logical index array..
LIArrayz=cell(length(Eqlz),1);

for ii =1:length(Eqlz)

    EqTmp=Eqlz{ii}; %save equality statement string as "EqTmp" for this pass
    
    % Iterate through the names in itParNamz checking whether each is
    % present in EqTmp, and if so replace it with the expression
    % referring to it's corresponding row in the itrParsCom matrix..
    for vv = 1:length(itParNamz)
        EqTmp = strrep(EqTmp,itParNamz{vv},strcat('itrParsCom(',num2str(vv),",:)"));
    end
    
    % Evaluate the updated EqTmp statement as a command: "idx = EqTmp" to 
    % generate a logical index array corresponding to columns in itrParsCom 
    % meeting the equality statement..
    
    %Build Command
    cmd=strcat("idx = ",EqTmp,";"); 
    
    % Replace multiplication and division symbols (if present) with dot
    % notation equivilants for matrix operation..
    cmd = strrep(cmd,"*",".*");
    cmd = strrep(cmd,"/","./");
    
    eval(cmd) %now run it..

    % save logical index array output "idx" to LIArrayz...
    LIArrayz(ii,1)={idx};

end
                               
% Now multiply the logical index arrays in LIArrayz to combine them..
% Multiply all the logical index arrays using the "all" function
idxCombined = all(cell2mat(LIArrayz), 1);

% Finally, extract the columns which met all the conditions in Eqlz
% logically indexing itrParsCom with idxCombined..
itrParsCom_matches = itrParsCom(:, idxCombined);

end