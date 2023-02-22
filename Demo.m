% This code is a companion to our submitted ICIP 2023 paper for algorithms
% APIRL1_AM and PIRL1_AM

clc
clear all;
close all; 

imageName = 'pepper512.png';

Img = double(imread(imageName)); %Your Image goes here

N = numel(Img);

[row, col] = size(Img);

row = int2str(row);
col = int2str(col);

imageSize = [row 'x' col];

K = fspecial('gaussian', [17 17], 7); 
f = imfilter(Img,K,'circular');


BSNR = 30;
sigma = BSNR2WGNsigma(f, BSNR);

fprintf('The noise std of the observed image: %g.\n', sigma); 


f = f +  sigma * randn(size(Img)); %Add a little noise

%**************Initialize parameters for denoising*****************

opts.mu           = 30;  %smaller the more noise filtering
opts.beta         = 0.009 ; % for small BSNR, use larger beta values 
opts.rho          = 1;
opts.Nit           = 1000;
opts.tol           = 1.0e-8;
opts.p            = 0.1;

%%% For Accelerated non-convex AM %%%%
optsNcvxAcc.mu = 30;
optsNcvxAcc.beta = 0.009;


%%%% Test algorithms %%%%%

out = APIRL1_AM(f, Img, K, opts, optsNcvxAcc);
out2 =PIRL1_AM(f, Img, K, opts);


figure;
imshow(out.sol,[]);
title(sprintf('APIRL1-AM (PSNR = %3.3f dB,SSIM = %3.3f, cputime %.3f s) ',...
                       psnr_fun(out.sol,Img),ssim_index(out.sol,Img), out.cpuTime));


figure;
imshow(out2.sol,[]);
title(sprintf('PIRL1-AM (PSNR = %3.3f dB,SSIM = %3.3f, cputime %.3f s) ',...
                       psnr_fun(out2.sol,Img),ssim_index(out2.sol,Img), out2.cpuTime));
                 
                   
figure;
imshow(uint8(f));
title(sprintf('Noisy (PSNR = %3.3f dB, SSIM = %3.3f)',psnr_fun(f,Img), ssim_index(f,Img)));

figure;
imshow(Img,[]);
title('Original');

figure;
 
semilogy(out.relativeError,'Linewidth',2.5,'Color','black');hold
axis tight;

semilogy(out2.relativeError,'--','Linewidth',3.0,'Color','red');

xlabel('Iterations (k)','FontSize',25,'interpreter','latex');
ylabel('Relative Error','FontSize',25,'interpreter','latex');
axis tight;
grid on;

l = legend('APIRL1-AM','PIRL1-AM','Acc AM');
set(l,'interpreter','latex','FontSize', 25);
set(gca, 'FontSize',20)
