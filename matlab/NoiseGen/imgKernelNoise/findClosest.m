function [closestVal,idx] = findClosest(target, array)
    % Calculate the absolute difference between target and each value in array
    diffs = abs(array - target);
    
    % Find the index of the smallest difference
    [~, idx] = min(diffs);
    
    % Return the closest value
    closestVal = array(idx);
end
