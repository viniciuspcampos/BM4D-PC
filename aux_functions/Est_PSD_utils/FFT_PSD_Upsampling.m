function glob_FFT_PSD = FFT_PSD_Upsampling(loc_FFT_PSD, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  FFT_PSD_Upsampling converts the local FFT-PSD ("loc_FFT_PSD") into
%  a global FFT-PSD ("glob_FFT_PSD") of size "N" via upsampling.
%
%  The upsampling is based on solving a least-squares optimization, minimizing
%  the error in the local FFT-PSD given its global form, while enforcing non-negativity,
%  smoothness, and symmetries of global Fourier-domain PSDs. We solve the optimization 
%  iteratively as in Section 6.3.2 of:
%           L. Azzari, L. R. Borges, and A. Foi, "Modeling and estimation of signal  
%           dependent and correlated noise," in Denoising of Photographic Images
%           and Video. Springer, 2018, pp. 1–36.
%  
%  
%
%  FUNCTION INTERFACE:
%     glob_FFT_PSD = FFT_PSD_Upsampling(loc_FFT_PSD, N)
% 
%  INPUTS:
%     1) loc_FFT_PSD          (2D/3D array) : Local FFT-PSD
%     2) N                    (double)      : the size of output PSD
%
%  OUTPUTS:
%     1) glob_FFT_PSD         (2D/3D array) : Global FFT-PSD (upsampled PSD)
%
%
%  EXAMPLES:
%  prerequisite : Est_rootPSD.m (to compute the local FFT-PSD using the sample variance)
%%-------------------------------------------------------------------------
% g1                   = fspecial('gaussian',7,1.5); % a lowpass Gaussian filter 
% g2                   = fspecial('gaussian',7,0.5); % a lowpass Gaussian filter
% kernel1              = g1-g2;                      % a bandpass Gaussian filter 
%%-------------------------------------------------------------------------
% kernel2              = rand(13);                   % a random kernel
%%-------------------------------------------------------------------------
% kernel3              = .1*randn(size(kernel1))+kernel1; % another kernel
%%-------------------------------------------------------------------------
%
% g                    = kernel3; % kernel1 % kernel2; % selected correlation kernel
% 
% N                    = [1023 511];                 % N1xN2: size of 2-D noise signal
% L                    = 11;                         % LxL: size of converted 2-D PSD
% w                    = randn(N);                   % a standard white Gaussian noise of size N
% eta                  = imfilter(w,g,'circular');   % a correlated noise
% 
% figure, imagesc(eta), colormap('gray'), axis off
% 
% glob_FFT_PSD_kernel   = ...
%             abs(fft2(g,N(1),N(2))).^2*(prod(N));   % global FFT-PSD of "eta"
%                                                    % computed from the correlation kernel
% loc_FFT_PSD       = (Est_rootPSD(eta, [L L],...
%              {'fft','fft'}, 0,'sVar')).^2;         % local FFT-PSD of "eta"
%                                                    % computed using sample variance
% 
% glob_FFT_PSD = FFT_PSD_Upsampling(loc_FFT_PSD,N);
%                                                    % estimated global FFT-PSD of "eta"
%                                                    % via upsampling
% 
% figure, 
% subplot(1,3,1), imagesc(fftshift(glob_FFT_PSD_kernel)), colormap('jet'),
% colorbar, title('the global FFT-PSD computed from kernel')
%         subplot(1,3,2), imagesc(fftshift(glob_FFT_PSD)), colormap('jet'),
% colorbar, title('the global FFT-PSD via upsampling')
%         subplot(1,3,3), imagesc(fftshift(loc_FFT_PSD)), colormap('jet'),
% colorbar, title('the local FFT-PSD')
%--------------------------------------------------------------------------
%
%               AUTHORS:  Alessandro Foi & Nasser Eslahi
%                 contact: firstname.lastname@tuni.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2006-2022 Tampere University.
% All rights reserved.
% This work (software, material, and documentation) shall only be used for nonprofit 
% noncommercial purposes.
% Any unauthorized use of this work for commercial or for-profit purposes is prohibited.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N1 = N(1);
N2 = N(2);
if numel(N)==2
    N3 = 1;
elseif numel(N)==3
    N3 = N(3);
else
    error('this version only accepts 2D and 3D PSDs!')
end


n1         = size(loc_FFT_PSD,1);
Basis1     = zeros(N1,n1);
for j1=1:n1
    eee=zeros(n1,1); eee(j1)=1;
    Basis1(1:n1,j1)=ifft(eee); 
    
end
Phi2_1     = abs(fft(Basis1)).^2/N1^2;       %% frame
Phi2dual_1 = Phi2_1*pinv(Phi2_1'*Phi2_1);    %% dual frame
clear Basis1


n2         = size(loc_FFT_PSD,2);
Basis2     = zeros(N2,n2);
for j2=1:n2
    eee=zeros(n2,1); eee(j2)=1;
    Basis2(1:n2,j2)=ifft(eee);  
end
Phi2_2     = abs(fft(Basis2)).^2/N2^2;       %% frame
Phi2dual_2 = Phi2_2*pinv(Phi2_2'*Phi2_2);    %% dual frame
clear Basis2



n3         = size(loc_FFT_PSD,3);
Basis3     = zeros(N3,n3);
for j3=1:n3
    eee=zeros(n3,1); eee(j3)=1;
    Basis3(1:n3,j3)=ifft(eee);  
end
Phi1_3     = abs(fft(Basis3)).^2/N3^2;       %% frame
Phi1dual_3 = Phi1_3*pinv(Phi1_3'*Phi1_3);    %% dual frame
clear Basis3


glob_FFT_PSD  = zeros(N1,N2,N3);
glob_residual = glob_FFT_PSD;
loc_FFT_PSD   = reshape(loc_FFT_PSD,[n1*n2 n3]);

Niter = 100;

if N3==1
    % this is faster in 2D cases 
    pad = @(I,d)[I(end-d(1)+1:end,end-d(2)+1:end,:)   I(end-d(1)+1:end,:,:)   I(end-d(1)+1:end,1:d(2),:)
                 I(:,end-d(2)+1:end,:)                I                       I(:,1:d(2),:)
                 I(1:d(1),end-d(2)+1:end,:)           I(1:d(1),:,:)           I(1:d(1),1:d(2),:)];
else
    % this is faster in 3D cases 
    % (if image processing toolbox is not installed, you can use the above "pad" handle_function)
    pad = @(I,d) padarray(I, d, 'both','circular'); 
end

smth_krnl=ones(5); smth_krnl=smth_krnl/sum(smth_krnl(:)); % smoothing kernel

for jj=1:Niter
        
    rho = 1;  %% step-size
    
    dummy            = zeros(n1,n2,N3);
    %%% projection onto desirable subspace 
    for i=1:N3
        dummy(:,:,i) = Phi2_1'*glob_FFT_PSD(:,:,i)*Phi2_2;
    end
    dummy            = reshape(dummy,[n1*n2 N3]);
    est_loc_FFT_PSD  = dummy*Phi1_3;

    
    loc_residual     = loc_FFT_PSD - est_loc_FFT_PSD;     %% local FFT-PSD residual
    
    %%% backward projection by dual frames
    cnvt_1d_residual = reshape( loc_residual*Phi1dual_3', [n1 n2 N3]);
    for i=1:N3
        glob_residual(:,:,i) = Phi2dual_1*cnvt_1d_residual(:,:,i)*Phi2dual_2'; %% global FFT-PSD residual
    end
    
    glob_FFT_PSD     = glob_FFT_PSD + rho* glob_residual; %% adding global residual to the previous estimate
    
    %%% constraints
  %  glob_FFT_PSD(2:end,1,:)     = (glob_FFT_PSD(2:end,1,:)+glob_FFT_PSD(end:-1:2,1,:))/2; %% impose symmetry
  %  glob_FFT_PSD(1,2:end,:)     = (glob_FFT_PSD(1,2:end,:)+glob_FFT_PSD(1,end:-1:2,:))/2; %% impose symmetry
  %  glob_FFT_PSD(2:end,2:end,:) = (glob_FFT_PSD(2:end,2:end,:)+glob_FFT_PSD(end:-1:2,end:-1:2,:))/2; %% impose symmetry

    
    glob_FFT_PSD = max(0,glob_FFT_PSD); %% impose non-negativity (it's a PSD!)
    
    
    for juu=1:2  
        glob_FFT_PSD = pad(glob_FFT_PSD,(size(smth_krnl)-1)/2);
        glob_FFT_PSD = convn(glob_FFT_PSD, smth_krnl,'valid'); %% impose smoothness      

%%% glob_FFT_PSD = imfilter(glob_FFT_PSD, smth_krnl, 'circular'); % equivalent but slower

    end
 
    
end

glob_FFT_PSD = glob_FFT_PSD./numel(loc_FFT_PSD)^2; 