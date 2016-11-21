%% Image segmentation and extraction
%% Read Image
width = 400;
imagen=imread('test2.jpg');
imagen = imresize(imagen,[width,width]);



%% Show image
figure(1)
imshow(imagen);
title('INPUT IMAGE WITH NOISE')


% Get the black white version
imagen = im2bw(imagen,0.75);
imshow(imagen)



%% Remove all object containing fewer than 30 pixels
imagen = bwareaopen(imagen,30);
pause(1)
%% Show image binary image
figure(2)
imshow(~imagen);
title('INPUT IMAGE WITHOUT NOISE')
%% Label connected components
[L Ne]=bwlabel(imagen);
%% Measure properties of image regions
propied=regionprops(L,'BoundingBox');
hold on
%% Plot Bounding Box
for n=1:size(propied,1)
    rectangle('Position',propied(n).BoundingBox,'EdgeColor','g','LineWidth',2)
end
hold off
pause (1)
%% Objects extraction
figure
for n=1:Ne
    [r,c] = find(L==n);
    n1=imagen(min(r):max(r),min(c):max(c));
    imshow(~n1);
    pause(0.5)
end