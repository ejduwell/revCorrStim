function [firstPart, secondPart] = splitTitle(stringIn, L, splitChar)
    % Find the indices of the splitChar in the string
    splitIndices = strfind(stringIn, splitChar);
    
    if isempty(splitIndices)
        error('The specified splitChar was not found in the input string.');
    end
    
    % Find the index of the nearest splitChar to L
    [~, nearestIdx] = min(abs(splitIndices - L));
    splitIdx = splitIndices(nearestIdx);
    
    % Split the string into two parts
    firstPart = stringIn(1:splitIdx - 1);
    secondPart = stringIn(splitIdx + 1:end);
end