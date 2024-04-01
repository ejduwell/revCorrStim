function nearestPowerOfTwo = findNearestPowerOfTwo(number)
    if number <= 0
        error('Number must be greater than zero.');
    end
    
    nearestExponent = round(log2(number));
    nearestPowerOfTwo = 2^nearestExponent;
end
