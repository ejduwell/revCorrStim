function errorSignal= compareFFTs(im1, im2)

    %im1FFT = fft2(im1);
    %im2FFT = fft2(im2);
%     im1FFT_frq = real(im1FFT);
%     im1FFT_phz = imag(im1FFT);
%     im2FFT_frq = real(im2FFT);
%     im2FFT_phz = imag(im2FFT);
 


    % Compute power spectrum/fft for image 1
    F1 = fft2(double(im1));
    F1 = fftshift(F1);
    S1 = abs(F1).^2;
    %max1 = max(S1, [], 'all');
    
    % Compute the power spectrum/fft for image 2
    F2 = fft2(double(im2));
    F2 = fftshift(F2);
    S2 = abs(F2).^2;
    %max2 = max(S2, [], 'all');
    
    %max12 = max([max1,max2]);

    % normalize S1 and S2 based on maximum value present in either of the 
    % two spectra..
%     S1 = S1./max12;
%     S2 = S2./max12;
      %errorSignal= mean2((S1-S2).^2);

    SpercDiff = abs(S2-S1)./(S1+S2);
    errorSignal = mean2(SpercDiff.^2);
    if isnan(errorSignal)
        SpercDiff(isnan(SpercDiff))=0;
        errorSignal = mean2(SpercDiff.^2);
    end
    %errorSignal= mean2((S1-S2).^2);

    
%     errorSignal_frq = (mean2((im1FFT_frq-im2FFT_frq).^2));
%     errorSignal_phz = (mean2((im1FFT_phz-im2FFT_phz).^2));
    %pause = "";
end