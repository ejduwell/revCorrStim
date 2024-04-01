function result = isPowerOfTwo(number)
    result = number > 0 && bitand(number, number - 1) == 0;
end
