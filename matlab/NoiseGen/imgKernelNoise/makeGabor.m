function [MyGab,Gaussian,mymask] = makeGabor(ori, c, nPix, ar, sf, sigma, wDeg, expMask)


useDefault = 0;
rsData = []; %[-1 1];

if useDefault
    c           = 1;  % Contrast
    ar          = 1;   % Aspect ration
    wDeg        = 1;   % size of image (in degrees)
    nPix        = 200; % resolution of image (pixels);
    sigma       = .15; % width of Gaussian (1/e half-width)
    ori         = 90;  % deg (counter-clockwise from horizontal)
    sf          = 4;   % spatial frequency (cycles/deg)
end

[x,y] = meshgrid(linspace(-wDeg/2,wDeg/2,nPix+1));
x = x(1:end-1,1:end-1);
y = y(1:end-1,1:end-1);


Gaussian = exp(-(x.^2/(sigma*ar)^2 + y.^2/(sigma*ar)^2));
if ~isempty(rsData)

Gaussian = rescale(Gaussian,rsData(1),rsData(2));

end

ramp     = cos(ori*pi/180)*x - sin(ori*pi/180)*y;
grating  = c*sin(2*pi*sf*ramp);
MyGab    = grating.*Gaussian;
mymask   = Gaussian;

if expMask ==1
    
    a = figure();
    [C,~] = contour(mymask);
    close(a);
    breaks = find(C(1,:)<=1);
    C = round(C);
    xx = C(1,2:breaks(2)-1);
    yy = C(2,2:breaks(2)-1);
    mymask = zeros(nPix);
    for xdim = 1:nPix
        for ydim = 1:nPix
            [in,on] = inpolygon(xdim,ydim,xx,yy);

            if in || on
                mymask(xdim,ydim) = 1;
            end
        end
    end
end

end

