%%%read the result
denoised_image = csvread('../data/temp.csv');
noise_image    = csvread('../data/flower_256_noise.csv');

h = size(denoised_image,1);

%%%delete the last column which matlab adds and its value is zero
%%%everywhere

denoised_image = denoised_image(1:h,1:h);

%plot image
figure(1);
imshow(denoised_image);
title('denoised image');

%plot residue
figure(2);
imagesc(denoised_image-noise_image); 
colormap gray;
title('residue');

figure(3);
imshow(noise_image);
title('noisy image');