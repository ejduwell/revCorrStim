function coefficient = myBandPassFilter(distance,lowFreq,highFreq)
%     lowFreq = 30; % Low frequency cutoff
%     highFreq = 60; % High frequency cutoff
% 
%     lowFreq = 0; % Low frequency cutoff
%     highFreq = 10; % High frequency cutoff

    if distance >= lowFreq && distance <= highFreq
        coefficient = 1;
    else
        coefficient = 0;
    end
end
