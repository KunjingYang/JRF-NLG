clear
addpath('JRFNLG_toolbox')
true_image =load('.\data\Superballs.mat');
true_image  = true_image.S ;
true_image =double(true_image ); 
true_image =true_image /max(true_image(:));
par.true_image = true_image;
[M,N,L]=size(true_image );


sf =16;  s0=1;  sz=[M N];
psf        =    fspecial('gaussian',7,2);
par.fft_B      =    psf2otf(psf,sz);
par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);

F = create_F();

 for band = 1:size(F,1)
        div = sum(F(band,:));
        for i = 1:size(F,2)
            F(band,i) = F(band,i)/div;
        end
 end
 
S_bar = hyperConvert2D(true_image );
hyper= par.H(S_bar); % Ä£ºý+½µ²ÉÑù ¼´BS
MSI = hyperConvert3D((F*S_bar), M, N);
HSI =hyperConvert3D(hyper,M/sf, N/sf );

% SNRh=35; 
% hyper= par.H(S_bar);
% sigmah =sqrt(sum(hyper(:).^2)/(10^(SNRh/10))/numel(hyper));
% rng(10,'twister')
% hyper= hyper+sigmah*randn(size(hyper));
% HSI=hyperConvert3D(hyper,M/sf,N/sf);
% 
% Y=F*S_bar;
% SNRm=40; 
% sigmah1 =sqrt(sum(Y (:).^2)/(10^(SNRm/10))/numel(Y ));
% rng(10,'twister')
% MSI = hyperConvert3D((Y+sigmah1*randn(size(Y))), M, N);
 

nlrgs.K = 210; nlrgs.patchsize = 7; nlrgs.p1 = 11; nlrgs.p2 = 5; nlrgs.tao = 5; nlrgs.mu = 1e-5;
nlrgs.alpha=6e-3;  nlrgs.beta = 1e-5; nlrgs.mul = 2e-3; nlrgs.mue = 5e-3; nlrgs.theta = 9;

rafpara.iscal = 1;
t0=clock;
[X,~,kappa] = NLRGS_Fus(HSI,MSI,F,sf,nlrgs,par,5,16,rafpara);
t4=etime(clock,t0);
[PSNR,RMSE,ERGAS,SAM, ~,SSIM,~,~] = quality_assessment(double(im2uint8(true_image)),double(im2uint8(X)),0,1/sf);


