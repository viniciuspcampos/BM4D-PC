
function [I_denoised]= BM4D_PC(I_Noisy, Estimated_NoiseMap, glob_PSD, profile)


if(~exist('profile', 'var'))
    profile = 'np';
end

%****** Deviding noisy image by the Estimated map ****%
[X, Y, Z, W] = size(I_Noisy);

I_Noisy = I_Noisy./Estimated_NoiseMap;


% ****** PCA decomposition for denoising **********%


I_Noisy = reshape(I_Noisy,[X*Y*Z W]);

I_noisy_cov = I_Noisy'*I_Noisy;


[~, ~, v] = svd(I_noisy_cov, 'econ');


I_Noisy = I_Noisy*v;

I_Noisy = reshape(I_Noisy,[X Y Z W]);


%% ** adjust psd

glob_PSD = glob_PSD /  (norm(glob_PSD(:),1)/numel(glob_PSD)^2);  % Just enforcing/guaranteeing it is normalized (var = 1)


variances = ones(1,W); % due to the normalization by Estimated map, variance is one on every PC


%% bm4d profile setup
BM4D_profile= get_BM4D_profile(profile); %get profile settings

stage_arg = BM4DProfile.ALL_STAGES; %perform both Hard-Thresholding and Wiener stages of BM4D


%% Denoise first Component and getting block-matching

fprintf('Denoising first Component and getting block-matching poisitions...\n')

I_denoised = double(I_Noisy);

i=1;

tic_first_pc = tic;

[I_denoised(:,:,:,i), match_arrs] = BM4D(I_Noisy(:,:,:,i),(glob_PSD.*variances(i)), BM4D_profile, stage_arg, {true, true});

elap_time_first_pc = toc(tic_first_pc);

fprintf('Elapsed time denoising first PC = %.2f (sec) ...\n',elap_time_first_pc)

%% Denoise remaining PCs
fprintf('Denoising remaining Components ...\n')
fprintf('Estimated time to finish = %.2f (sec) ...\n',elap_time_first_pc*(W-1));

for i=2:W

    [I_denoised(:,:,:,i)] = BM4D(I_Noisy(:,:,:,i), (glob_PSD.*variances(i)), BM4D_profile, stage_arg, match_arrs);
            
end


%%  'Inverse' PCA

I_denoised = reshape(I_denoised,[X*Y*Z W]);

I_denoised = I_denoised*v' ;

I_denoised = reshape(I_denoised,[X Y Z W]);


%******* Multiplying denoised image by the Estimated map ****%
I_denoised = I_denoised.*Estimated_NoiseMap;


end
