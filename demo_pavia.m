clear; %clc ;

% addpath
addpath JRFNLG_toolbox;
addpath data ;
addpath results ;


%% Import data to obtain the initial HSI and MSI
true_image   =   imread('.\data\original_rosis.tif');
true_image   =   true_image(1:256,1:256,11:end);


true_image   =   double(true_image); 
true_image   =   true_image/max(true_image(:));
par.true_image = true_image;

[M,N,~]      =   size(true_image);
sf = 8;   s0 = 1;   sz = [M N];   par.sf = sf;

psf            =    fspecial('gaussian',7,2);          % blur kernel
par.fft_B      =    psf2otf(psf,sz);
par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 ); % Spatial downsampling operator
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);

load('.\data\R.mat');
F=R; F = F(:,6:end-5);                          
for band = 1:size(F,1)
        div = sum(F(band,:));
        for i = 1:size(F,2)
            F(band,i) = F(band,i)/div;
        end
end
par.R = F;                                             % Spectral dowmsamping operator

% obtain the registered HSI and MSI, and add noise
S_bar = hyperConvert2D(true_image);
SNRh=35;
hyper= par.H(S_bar);
sigmah =sqrt(sum(hyper(:).^2)/(10^(SNRh/10))/numel(hyper));
rng(10,'twister')
hyper= hyper+sigmah*randn(size(hyper));
HSI=hyperConvert3D(hyper,M/sf,N/sf);

Y=F*S_bar;
SNRm=40; 
sigmah1 =sqrt(sum(Y (:).^2)/(10^(SNRm/10))/numel(Y ));
rng(10,'twister')
MSI = hyperConvert3D((Y+sigmah1*randn(size(Y))), M, N);


%% Transform HSI to generate unregistered HSI
unregisHsi = HSI;

unregisHsi = imrotate(unregisHsi,10);                 % rotation
% unregisHsi = imtranslate(unregisHsi,[5, 5]);          % translation
% unregisHsi = imresize(unregisHsi,1.1);                  % scaling
% unregisHsi = barrel_distortion (unregisHsi, 1e-3);    % barrel_distortion  

% tran_matrix=[1 0.5 0; 0.5 1 0; 0 0 1];                % flip
% tform=maketform("affine",tran_matrix);
% unregisHsi=imtransform(unregisHsi,tform);

% Restore to original size
unregisHsi = unregisHsi(1:M/sf,1:M/sf,:); 

%% Save as batch images
Spec_bands = size(unregisHsi, 3);
str = '.\data\HSI\';
points = [4.1500, 45.6800; 19.6000, 19.6000];
for band = 1:Spec_bands   % mat to bmp
    if band < 10
        name = sprintf('00%d.bmp',band);
        imwrite(unregisHsi(:,:,band),[str,name] );
        name_mat = sprintf('00%d-points.mat',band);
        filename = fullfile(str, name_mat);  
        save(filename, 'points');   
    else
        name = sprintf('0%d.bmp',band);
        imwrite(unregisHsi(:,:,band),[str,name] );
        name_mat = sprintf('0%d-points.mat',band);
        filename = fullfile(str, name_mat);  
        save(filename, 'points'); 
    end
end
par.unregisHsi = unregisHsi;

%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data') ;
pointPath = fullfile(currentPath,'data') ; % path to files containing initial feature coordinates
userName = 'HSI' ;

% output path
destRoot = fullfile(currentPath,'results') ;
destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

%% define parameters

% display flag
rafpara.DISPLAY = 0 ;

% save flag
rafpara.saveStart = 1 ;
rafpara.saveEnd = 1 ;
rafpara.saveIntermedia = 0 ;

rafpara.canonicalImageSize = [ M/sf,  M/sf  ];     
rafpara.canonicalCoords = [ 0.5   44 ; 17 17  ]; 

% parametric tranformation model
rafpara.transformType = 'SIMILARITY'; 
% one of 'TRANSLATION', 'SIMILARITY', 'AFFINE','HOMOGRAPHY'

rafpara.numScales = 1 ; % if numScales > 1, we use multiscales


%% Get training images
% get initial transformation
transformationInit = 'SIMILARITY'; 
[fileNames, transformations, numImages] = get_training_images( imagePath, pointPath, userName, rafpara.canonicalCoords, transformationInit, Spec_bands) ;

%% JRF: Joint registration and fusion
rafpara.step = 1;     % Iteration number of JRF model
rafpara.subspaceDim  = 3;  % subspace dimension
rafpara.tao = 10;   
rafpara.mu1 = 10;     
rafpara.mu2 = 10;     
rafpara.iscal = 1;    % Whether calculate PSNR in each iteration

[D, A, xi,t1, Z, kappa ] = JRF_outer(fileNames, transformations, numImages, rafpara, destDir, MSI, par);

%% NLG: Nonlocal low-rank and group sparse
load('Registered_HSI.mat');  Registered_HSI = MN_3D;
nlrgs.K=200;  nlrgs.patchsize=7;  nlrgs.p1=3;  nlrgs.p2=90; nlrgs.tao=15; nlrgs.mu=1e-5; m=1;
nlrgs.alpha=2e-3; nlrgs.beta=2e-3; nlrgs.mul=5e-3; nlrgs.mue=5e-3; nlrgs.theta=7; n=3;

[X, t2] = NLRGS_Fus(Registered_HSI, MSI, par.R, sf, nlrgs, par, m,n,rafpara);

%% Result of RAF-NLRGS

[PSNR,RMSE,ERGAS,SAM, ~,SSIM,~,~] = quality_assessment(double(im2uint8(par.true_image)),double(im2uint8(X)),0,1/sf);
the_index_of_RAF_NLRGS = ["PSNR" "SSIM" "ERGAS" "SAM"  "RMSE" "time"; PSNR SSIM ERGAS SAM RMSE t1+t2]

subplot(131);   imshow(func_hyperImshow(unregisHsi,[10,20,30]));             title('Initial LR-HSI')
subplot(132);   imshow(func_hyperImshow(Registered_HSI,[10,20,30]));         title('Registrated LR-HSI')
subplot(133);   imshow(func_hyperImshow(X,[10,20,30]));                      title('Fused HR-HSI')

