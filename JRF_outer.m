function [D, A, xi, timeConsumed,Z,kappa1 ] = JRF_outer(fileNames, transformations, ...
                                                        numImages, rafpara, destDir, MSI, par)
disp(' JRF model is running...' );
%% read and store full images

fixGammaType = 1 ;

if ~fixGammaType
    if exist(fullfile(rootPath, userName, 'gamma_is_ntsc'), 'file')
        gammaType = 'ntsc' ;
    elseif exist(fullfile(rootPath, userName, 'gamma_is_srgb'), 'file')
        gammaType = 'srgb' ;
    elseif exist(fullfile(rootPath, userName, 'gamma_is_linear'), 'file')
        gammaType = 'linear' ;
    else
        error('Gamma type not specified for training database!  Please create a file of the form gamma_is_*') ;
    end
else
    gammaType = 'linear' ;
end

sigma0 = 1/5 ;
sigmaS = 1 ;

deGammaTraining = true ;

I0 = cell(rafpara.numScales,numImages) ; % images
I0x = cell(rafpara.numScales,numImages) ; % image derivatives
I0y = cell(rafpara.numScales,numImages) ;

for fileIndex = 1 : numImages
    
    currentImage = double(imread(fileNames{fileIndex}));
    
    % Use only the green channel in case of color images    
    if size(currentImage,3) > 1,   currentImage = currentImage(:,:,2);            end
    if deGammaTraining,      currentImage = gamma_decompress(currentImage, gammaType); end
    
    currentImagePyramid = gauss_pyramid( currentImage, rafpara.numScales,...
        sqrt(det(transformations{fileIndex}(1:2,1:2)))*sigma0, sigmaS );
        
    for scaleIndex = rafpara.numScales:-1:1
        I0{scaleIndex,fileIndex} = currentImagePyramid{scaleIndex};
        
        % image derivatives
        I0_smooth = I0{scaleIndex,fileIndex};
        I0x{scaleIndex,fileIndex} = imfilter( I0_smooth, (-fspecial('sobel')') / 8 );
        I0y{scaleIndex,fileIndex} = imfilter( I0_smooth,  -fspecial('sobel')   / 8 );
    end   
end



%% get the initial input images in canonical frame

imgSize = rafpara.canonicalImageSize ; 

xi_initial = cell(1,numImages) ; % initial transformation parameters
for i = 1 : numImages
    if size(transformations{i},1) < 3
        transformations{i} = [transformations{i} ; 0 0 1] ;
    end
    xi_initial{i} = projective_matrix_to_parameters(rafpara.transformType,transformations{i});
end

%D = [] ;
D= zeros(imgSize(1)*imgSize(2), numImages);
for fileIndex = 1 : numImages
    % transformed image
    Tfm = fliptform(maketform('projective',transformations{fileIndex}'));
            
    I   = vec(imtransform(I0{1,fileIndex}, Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
    y   = I; 
    y = y / norm(y) ; % normalize

    %D = [D y] ;
    D(:,fileIndex) = y ;
end

if rafpara.saveStart
    save(fullfile(destDir, 'original.mat'),'D','xi_initial') ;
end


%% start the main loop

frOrig = cell(1,numImages) ;
T_in = cell(1,numImages) ;

T_ds = [ 0.5,   0, -0.5; ...
         0,   0.5, -0.5   ];
T_ds_hom = [ T_ds; [ 0 0 1 ]];

numIterOuter = 0 ; 
numIterInner = 0 ;

tic % time counting start

for scaleIndex = rafpara.numScales:-1:1 % multiscale
    
    iterNum = 0 ;  % iteration number of outer loop in each scale
    
    imgSize = rafpara.canonicalImageSize / 2^(scaleIndex-1) ;    
    xi = cell(1,numImages) ;
    
    for fileIndex = 1 : numImages
             
        if scaleIndex == rafpara.numScales
            T_in{fileIndex} = T_ds_hom^(scaleIndex-1)*transformations{fileIndex}*inv(T_ds_hom^(scaleIndex-1)) ;
        else
            T_in{fileIndex} = inv(T_ds_hom)*T_in{fileIndex}*T_ds_hom ;
        end
        
        % for display purposes
        if rafpara.DISPLAY > 0
            fr = [1 1          imgSize(2) imgSize(2) 1; ...
                  1 imgSize(1) imgSize(1) 1          1; ...
                  1 1          1          1          1 ];
            
            frOrig{fileIndex} = T_in{fileIndex} * fr;
        end
        
    end
    
    D= zeros(imgSize(1)*imgSize(2), numImages);
    
    

    kappa1 = zeros(1,rafpara.step);
    O1   = zeros(size(D')); 
    MSI3 = Unfold(MSI,size(MSI),3);
    O2   = zeros(size(MSI3));
    mu1         = rafpara.mu1;   
    mu2         = rafpara.mu2;
    
%% Outer loop
    for outer_iteration = 1:rafpara.step
        fprintf('Outer iteration: %d  ',outer_iteration);
        iterNum = iterNum + 1 ;
        numIterOuter = numIterOuter + 1 ;
        
        J = cell(1,numImages) ;
        
        for fileIndex = 1 : numImages

            % transformed image and derivatives with respect to affine parameters
            Tfm = fliptform(maketform('projective',T_in{fileIndex}'));
            
            % imtransform 实现图像的各种空间变换
            I   = vec(imtransform(I0{scaleIndex,fileIndex}, Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
            Iu  = vec(imtransform(I0x{scaleIndex,fileIndex},Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
            Iv  = vec(imtransform(I0y{scaleIndex,fileIndex},Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
            y   = I; %vec(I);

            Iu = (1/norm(y))*Iu - ( (y'*Iu)/(norm(y))^3 )*y ;
            Iv = (1/norm(y))*Iv - ( (y'*Iv)/(norm(y))^3 )*y ;

            %y = y / norm(y) ; % normalize
            % D = [D y] ;
            if outer_iteration==1
                D(:,fileIndex) = y ;
            end

            % transformation matrix to parameters
            xi{fileIndex} = projective_matrix_to_parameters(rafpara.transformType,T_in{fileIndex}) ; 
            
            % Compute Jacobian
            J{fileIndex} = image_Jaco(Iu, Iv, imgSize, rafpara.transformType, xi{fileIndex});
        end

        
        % JRF inner loop
        % -----------------------------------------------------------------
        % using QR to orthogonalize the Jacobian matrix
        for fileIndex = 1 : numImages
            [Q{fileIndex}, R{fileIndex}] = qr(J{fileIndex},0) ;
        end
        
        [D, A, delta_xi, numIterInnerEach,Z,AA] = JRF_inner(D, Q, MSI, par,O1, O2, rafpara);
        O1 = O1 + mu1*(par.H(AA)-D');
        O2 = O2 + mu2*(par.R*AA-MSI3);

        
        % stopping creterion
        ry  = par.R*D';
        zbs = par.H(Unfold(MSI,size(MSI),3));
        kappa1(outer_iteration) = norm(ry-zbs)/norm(zbs);
         
        
        for fileIndex = 1 : numImages
            delta_xi{fileIndex} = inv(R{fileIndex})*delta_xi{fileIndex} ;
        end
        numIterInner = numIterInner + numIterInnerEach ;

        % step in paramters
        for i = 1 : numImages           
            xi{i} = xi{i} +  delta_xi{i};
            T_in{i} = parameters_to_projective_matrix(rafpara.transformType,xi{i});
        end

    end
end

timeConsumed = toc;

disp(['time consumed by RAF is ' num2str(timeConsumed) ]);




