function [destinationRect] = scaleImgToScreen(Img, scale, centerPointX, centerPointY)

global EP

% Img is a matrix of img values
% Scale is an proportion (0-1) of the screen size

% Get the size of the image
[Y, X, ~] = size(Img);

% Get the aspect ratio of the image. We need this to maintain the aspect
% ratio of the image when we draw it different sizes. Otherwise, if we
% don't match the aspect ratio the image will appear warped / stretched
aspectRatio = X / Y;

% We will set the height of each drawn image 
% to a fraction of the screens height
imageHeight = EP.W_HEIGHT .* scale;
imageWidth = imageHeight .* aspectRatio;

% Make the destination rectangles for image
theRect = [0 0 imageWidth imageHeight];
destinationRect(:, 1) = CenterRectOnPointd(theRect, ...
    centerPointX, centerPointY);


end