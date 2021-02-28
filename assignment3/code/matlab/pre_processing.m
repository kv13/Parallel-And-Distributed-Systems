%%% read data
image_RGB = imread('../data/77_256.jpg');


%%% create greyscale image
I = rgb2gray(image_RGB);
imshow(I);
title('initial image');

%%% convert to double
d = im2double(I);

%%% add noise and save to csv format
noiseParams = {'gaussian', ...
                 0,...
                 0.001};

 

J = imnoise( d, noiseParams{:} );
csvwrite('../data/temp.csv',J);
fprintf("the image saved in folder data with name temp\n");