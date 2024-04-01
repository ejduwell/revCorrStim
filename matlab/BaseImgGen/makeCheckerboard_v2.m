function checkerboard = makeCheckerboard_v2(chkrBrdSz,sizeArrayOut)

% (ie this x2 will be the width/height of the array in units of total # of checks)
%chkrBrdSz=10; 
% square dimension of output array..
%sizeArrayOut=10000; 

% make initial chkrBrdSz*2 by chkrBrdSz*2 check array of 0s and 1s..
checkerboard = repmat(eye(2), chkrBrdSz, chkrBrdSz); 

% rescale to 0-255 and convert to uint8 format
checkerboard=uint8(rescale(checkerboard,0,255));

% resize the array to the desired dimensions
checkerboard=imresize(checkerboard,[sizeArrayOut,sizeArrayOut],"nearest");


end