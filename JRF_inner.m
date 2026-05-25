function [DD, A, dt_dual, iter,Z ,AA] = JRF_inner(D, J, MSI, par, lagM1, lagM2, rafpara)

subspaceDim = rafpara.subspaceDim;  
mu1         = rafpara.mu1;   
mu2         = rafpara.mu2;

alpha = 1e-3;   mu = 1e-4;   nu = 1e-2;  
tao   = 10;     theta = 8;   m = 5; n = 2;

R   = par.R;
sf  = par.sf;  
L   = size(D,2);
[nr, nc,l] = size(MSI); 

MSI3       = Unfold(MSI,size(MSI),3);
[dicSin,~] = svds(D', subspaceDim);  

e = eye(l);  F = zeros(size(MSI3));

HSI_ini  = hyperConvert3D(D', nr/sf, nc/sf); 
FBm      = par.fft_B;
HR_load1 = imresize(HSI_ini, sf,'bicubic'); % 初始点
FBmC = conj(FBm); FBCs = repmat(FBmC,[1 1 L]);
V2 = dicSin'*hyperConvert2D(HR_load1); A = V2;
n_dr = nr/sf; n_dc = nc/sf; O = zeros(size(V2));

dt_dual = cell(1,L) ;
dt_dual_matrix = zeros(nr*nc/(sf^2), L) ;

iter = 0;

FBs1  = repmat(FBm,[1 1 subspaceDim]);
FBCs1=repmat(FBmC,[1 1 subspaceDim]);


for k = 1:m
    iter = iter + 1;
    
    HSI3 = (D+dt_dual_matrix-lagM1'/mu1)';
    HSI = hyperConvert3D(HSI3, nr/sf, nc/sf); 
    HSI_int = zeros(nr,nc,L);
    HSI_int(1:sf:end,1:sf:end,:) = HSI; % 初始点
    HHH = ifft2((fft2(HSI_int).*FBCs));
    HSI_BS1 = hyperConvert2D(HHH);  % 即HSI_BS
    
    RdicSin = R*dicSin; DTD = dicSin'*dicSin; 
    CCC = DTD\(tao*RdicSin'*(e*MSI3+F-lagM2/mu2)+dicSin'*HSI_BS1); % （26）的第一项＋第三项
    C1 = DTD\(tao*(RdicSin'*RdicSin)+0.5*(nu+mu)*eye(size(dicSin,2)));  % H1
    [Q,Lambda] = eig(C1);
    Lambda = reshape(diag(Lambda),[1 1 subspaceDim]);
    InvLbd = 1./repmat(Lambda,[ sf*n_dr  sf*n_dc 1]);
    B2Sum = PPlus(abs(FBs1).^2./( sf^2),n_dr,n_dc);
    InvDI = 1./(B2Sum(1:n_dr,1:n_dc,:)+repmat(Lambda,[n_dr n_dc 1]));
    
    for i=1:n
        HR_HSI3 = 0.5*(nu*V2+mu*A+O); % （26）的第二项
        C3  = CCC+DTD\HR_HSI3;
        C30 = fft2(reshape((Q\C3)',[nr nc subspaceDim   ])).*InvLbd;
        temp = PPlus_s(C30/( sf^2).*FBs1,n_dr,n_dc);
        invQUF = C30-repmat(temp.*InvDI,[ sf  sf 1]).*FBCs1; 
        VXF    = Q*reshape(invQUF,[nr*nc subspaceDim])';
        A = reshape(real(ifft2(reshape(VXF',[nr nc subspaceDim ]))),[nr*nc subspaceDim])'; 
        
        B2 = A-O/(nu);
        % V2 = MCP_prox_tnn(B2,alpha/nu,theta,1);
        V2 =  Prox_CapL1(B2,alpha/nu,100);
%         V2=prox_tnn( B2, alpha/nu );
        %O=O+nu*(V2-A);
    end

    AA = dicSin*A;
    MN = par.H(AA);
    temp_T = D - MN';
    for i = 1 : L
        dt_dual{i} =  - J{i}'*temp_T(:,i) ; % delta_tau
        %dt_dual_matrix(:, i) = J{i}*dt_dual{i} ;
    end
    dt_dual_matrix = - temp_T;
    
%     e = (R*AA-F)/MSI3;
    
%     F = R*AA-e*MSI3;
end

%%
Z = hyperConvert3D(dicSin*A,nr, nc );
if rafpara.iscal == 1
    [PSNR,RMSE,ERGAS,SAM, ~,SSIM,~,~] = quality_assessment(double(im2uint8(par.true_image)),double(im2uint8(Z)),0,1/sf);
    fprintf('PSNR=%0.3f,SSIM=%0.3f, ERGAS=%0.3f, SAM=%0.3f, RMSE=%0.3f  \n',PSNR,SSIM,ERGAS,SAM,RMSE);
else
    fprintf('Here we do not calculate psnr. \n');
end


unregisHsi = par.unregisHsi;  
DD = D+dt_dual_matrix;
MN_3D = hyperConvert3D(DD', nr/sf, nc/sf );

save(['.\results\','Registered_HSI.mat'],'MN_3D');




