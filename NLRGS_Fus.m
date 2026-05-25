function  [Z,time, kappa, res] =  NLRGS_Fus( hsi, msi,R,sf,nlrgs,par,outer,inner,rafpara)

disp(' NLRGS model is running...' );

K=nlrgs.K; patchsize=nlrgs.patchsize; L1=nlrgs.p1; L2=nlrgs.p2; mu=nlrgs.mu;
alpha=nlrgs.alpha; beta=nlrgs.beta;  theta = nlrgs.theta;
mul = nlrgs.mul;  mue = nlrgs.mue; tao = nlrgs.tao;  overlap=4;

hsi3=Unfold(hsi,size(hsi),3);   msi3=Unfold(msi,size(msi),3);  L=size(hsi,3);
[lsm,~]=svds(hsi3,L1+L2);  dicl=lsm(:,1:L1);  dice = lsm(:,L1+1:L1+L2);  %C = dicl;

bparams.block_sz = [patchsize, patchsize];
bparams.overlap_sz=[overlap overlap];
[nr, nc,~]=size(msi); E = zeros(L2,nr*nc); F = E;

bparams.sz=[nr nc]; sz=[nr nc];
step=patchsize-overlap;
sz1=[1:step:sz(1)- bparams.block_sz(1)+1];
sz1=[sz1 sz(1)- bparams.block_sz(1)+1];
sz2=[1:step:sz(2)- bparams.block_sz(2)+1];
sz2=[sz2 sz(2)- bparams.block_sz(2)+1];
bparams.block_num(1)=length(sz1);
bparams.block_num(2)=length(sz2);

predenoised_blocks = ExtractBlocks1(msi, bparams);
Y2=Unfold(predenoised_blocks,size(predenoised_blocks),4);
aa=fkmeans(Y2,K);

hr_load1=imresize(hsi, sf,'bicubic'); % 初始点
FBm = par.fft_B; FBmC = conj(FBm); FBCs=repmat(FBmC,[1 1 L]);
V2=dicl'*hyperConvert2D(hr_load1); A = V2;  B =A;
n_dr=nr/sf; n_dc=nc/sf;  O=zeros(size(V2));

FBs1 = repmat(FBm,[1 1 L1]);  FBCs1 = repmat(FBmC,[1 1 L1]);
FBs2 = repmat(FBm,[1 1 L2]);  FBCs2 = repmat(FBmC,[1 1 L2]);
%%
tic
kappa = zeros(1,outer);
for k = 1:outer
    fprintf('Iteration: %d  ',k);
    hsi_int=zeros(nr,nc,L);
    MM = par.H(dice*E); Rdice = R*dice; 
    hsi_int(1:sf:end,1:sf:end,:) = hsi-hyperConvert3D(MM,n_dr, n_dc ); % 初始点
    HHH = ifft2((fft2(hsi_int).*FBCs));
    hsi_BS1 = hyperConvert2D(HHH);  % 即hsi_BS
    
    Rdicl = R*dicl; DTD = dicl'*dicl; 
    CCC = DTD\(tao*Rdicl'*(msi3-Rdice*E)+dicl'*hsi_BS1); % （26）的第一项＋第三项
    C1  = DTD\(tao*(Rdicl'*Rdicl)+0.5*(mul+mu)*eye(size(dicl,2)));  % H1
    [Q,Lambda] = eig(C1);
    Lambda = reshape(diag(Lambda),[1 1 L1]);
    InvLbd = 1./repmat(Lambda,[ sf*n_dr  sf*n_dc 1]);
    B2Sum  = PPlus(abs(FBs1).^2./( sf^2),n_dr,n_dc);
    InvDI  = 1./(B2Sum(1:n_dr,1:n_dc,:)+repmat(Lambda,[n_dr n_dc 1]));
    
    for i=1:inner
        hr_hsi3 = 0.5*(mul*V2+mu*A+O); % （26）的第二项
        C3   = CCC+DTD\hr_hsi3;
        C30  = fft2(reshape((Q\C3)',[nr nc L1   ])).*InvLbd;
        temp = PPlus_s(C30/( sf^2).*FBs1,n_dr,n_dc);
        invQUF = C30-repmat(temp.*InvDI,[ sf  sf 1]).*FBCs1; 
        VXF    = Q*reshape(invQUF,[nr*nc L1])';
        A = reshape(real(ifft2(reshape(VXF',[nr nc L1 ]))),[nr*nc L1])'; 
        
        B2 = A-O/(mul);
        %V2= MCP_prox_tnn(B2,alpha/mul,theta,1); 
        B2 = hyperConvert3D(B2,nr, nc );
        predenoised_blocks2 = ExtractBlocks1(B2, bparams);
        Z_block2=zeros( bparams.block_sz(1), bparams.block_sz(2),L1, bparams.block_num(1)* bparams.block_num(2));
        predenoised_blocks2=permute(predenoised_blocks2,[4 3 1 2]);
  
        for mn=1:max(aa)
            gg = find(aa==mn);
            XES = predenoised_blocks2(gg,:,:,:);
            [a, b, c, d ] = size(XES);
            XES = reshape(XES,[a b c*d]);
            % V22 = prox_tnn( XES, eta/2/mul );
            V22 = MCP_prox_tnn(XES,alpha/mul,theta,1); 
            V22 = reshape(V22,[a b c d]); 
            Z_block2(:,:,:,gg) = permute(V22,[3 4 2 1]);
        end  
        V2 = JointBlocks2(Z_block2, bparams);
        V2 =hyperConvert2D(V2);
        O  = O+mul*(V2-A);
    end
    kappa(k) = norm(A-B)/norm(B);  B = A;
    
    
    hsi_int = zeros(nr,nc,L);
    NN = par.H(dicl*A);
    hsi_int(1:sf:end,1:sf:end,:) = hsi-hyperConvert3D(NN,n_dr, n_dc ); 
    HHH = ifft2((fft2(hsi_int).*FBCs));
    hsi_BS2 = hyperConvert2D(HHH);  % 即hsi_BS
    DTD = dice'*dice; 
    CCC = DTD\(tao*Rdice'*(msi3-R*dicl*A)+dice'*hsi_BS2); % （26）的第一项＋第三项
    C1 = DTD\(tao*(Rdice'*Rdice)+0.5*(mue+mu)*eye(size(dice,2)));  % H1
    [Q,Lambda] = eig(C1);
    Lambda = reshape(diag(Lambda),[1 1 L2]);
    InvLbd = 1./repmat(Lambda,[ sf*n_dr  sf*n_dc 1]);
    B2Sum = PPlus(abs(FBs2).^2./( sf^2),n_dr,n_dc);
    InvDI = 1./(B2Sum(1:n_dr,1:n_dc,:)+repmat(Lambda,[n_dr n_dc 1]));
    for j = 1:10
        hr_hsi3 = 0.5*(mue*F+mu*E); 
        C3 = CCC+DTD\hr_hsi3;
        C30 = fft2(reshape((Q\C3)',[nr nc L2   ])).*InvLbd;
        temp  = PPlus_s(C30/( sf^2).*FBs2,n_dr,n_dc);
        invQUF = C30-repmat(temp.*InvDI,[ sf  sf 1]).*FBCs2;
        VXF    = Q*reshape(invQUF,[nr*nc L2])';
        E = reshape(real(ifft2(reshape(VXF',[nr nc L2]))),[nr*nc L2])'; 

        F = E;  vc = vecnorm(E);
        index = vc<sqrt(2*beta/mue);
        F(:,index) = 0;
    end

    Z=hyperConvert3D(dicl*A+dice*E, nr, nc );
    res=hyperConvert3D(dice*E, nr, nc );
if rafpara.iscal == 1
    [PSNR,RMSE,ERGAS,SAM, ~,SSIM,~,~] = quality_assessment(double(im2uint8(par.true_image)),double(im2uint8(Z)),0,1/sf);
    fprintf('PSNR=%0.3f,SSIM=%0.3f, ERGAS=%0.3f, SAM=%0.3f, RMSE=%0.3f  \n',PSNR,SSIM,ERGAS,SAM,RMSE);
else
    fprintf('Here we do not calculate psnr. \n');
end
end
time = toc;
disp(['time consumed by NLRGS is ' num2str(time) ]);

