function [] = ImageRecog()
% Given the image of a polynomial curve. Get the polynomial of degree n
% which best fits the curve.


width = 400;


% Get the image in RGB
RGB = imread('test2.jpg');
RGB = imresize(RGB,[width,width]);
imshow(RGB)

% Get the grayscale version
GS = rgb2gray(RGB);
imshow(GS)

% Get the black white version
BW = im2bw(GS,0.5);
imshow(BW)

% Invert the image Black White Complement
BWc = imcomplement(BW); 

% Find the white pixels
ii = find(BWc(:));

% Get the set of coordinate points of the white line.
[yy,xx]=ind2sub(size(BW),ii);

% Set the degree of the polynomial to fit the curve
n = 5;

% Fit the curve with a polynomial
p = polyfit(xx,yy,n);




% Use the polynomial to plot the curve
x1 = linspace(0,width,400);
y1 = polyval(p,x1);

figure
imshow(BW)
hold on
plot(x1,y1);
hold off


%dlmwrite('myFile.txt',[x1', y1'])
dlmwrite('myFile.txt',[xx,yy])




end