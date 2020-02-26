close all;
clear all;
clc;

% estimate reflectance and illumination
im=im2double(imread('75072.png'));
% You can modify the key parameters by intrinsicScene
[I,R]=intrinsicScene(im);

% show results
figure(1);
subplot(1,3,1);imshow(im);title('Input');
subplot(1,3,2);imshow(I);title('Illumination');
subplot(1,3,3);imshow(R,[]);title('Reflectance');
