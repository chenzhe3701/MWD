% chenzhe, 2018-03-18
% try FFT to filter out the high frequency part, see if grain boundaries
% can be better found.

close all;
addChenFunction;

I = imread('D:\p\m\DIC_Analysis\WE43_T6_C1_r0c0.tif');
% figure;
% imshow(I);

a = fftshift(fft2(I));
myplot(abs(a));

[nR,nC] = size(a);
rM = floor(nR/2)+1;
cM = floor(nC/2)+1;
sz = 100;
a(rM+(-sz:sz),cM+(-sz:sz)) = 0;
% myplot(abs(a));

I = ifft2(ifftshift(a));
figure;
imshow(I);