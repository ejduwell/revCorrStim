function pltGryImgsc(imgIn)
% Ethan got annoyed enough at imagesc for not plotting images with the
% proper dimensions/aspect ration that he wrote this.. It just plots an
% input image using imagesc() but sets colormap to "gray", ensures the
% aspect ratio is correct in the plot, and that a colorbar is included..

figure;
imagesc(imgIn);
colormap("gray");
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
colorbar;
%axis equal;

end