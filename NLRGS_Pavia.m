clear
load('.\data\R.mat');
addpath('JRFNLG_toolbox')

true_image=load('.\data\Pavia_Center.mat');
true_image=true_image.pavia;
true_image=true_image(1:256,1:256,10:end);
true_image=double(true_image); 
true_image=true_image/max(true_image(:));
par.true_image = true_image;
[M,N,L]=size(true_image);

sf =4;  s0=1;  sz=[M N];
psf        =    fspecial('gaussian',7,2);
par.fft_B      =    psf2otf(psf,sz);
par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);

F=R; F = F(:,6:end-5);
for band = 1:size(F,1)
        div = sum(F(band,:));
        for i = 1:size(F,2)
            F(band,i) = F(band,i)/div;
        end
end
S_bar = hyperConvert2D(true_image);
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


nlrgs.K=250; nlrgs.patchsize=6; nlrgs.p1=9; nlrgs.p2=40; nlrgs.tao=5; nlrgs.mu=1e-5; m=5;
nlrgs.alpha=1.7e-3; nlrgs.beta=2e-5; nlrgs.mul=1.4e-3; nlrgs.mue=2e-3; nlrgs.theta=8; n=15;



rafpara.iscal = 1;
t0=clock;
[X,~,kappa] = NLRGS_Fus(HSI,MSI,F, sf,nlrgs,par,m,n,rafpara);
t4=etime(clock,t0);


