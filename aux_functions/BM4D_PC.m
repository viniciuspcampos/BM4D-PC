
function [I_denoised]= BM4D_PC(I_Noisy, Estimated_NoiseMap, glob_PSD, profile)


if(~exist('profile', 'var'))
    profile = 'np';
end

%****** Deviding noisy image by the Estimated map ****%
[X, Y, Z, W] = size(I_Noisy);

I_Noisy = I_Noisy./Estimated_NoiseMap;


% ****** SVD for denoising **********%


I_Noisy = reshape(I_Noisy,[X*Y*Z W]);


[I_Noisy, vals, v] = svd(I_Noisy, 'econ');


I_Noisy = reshape(I_Noisy,[X Y Z W]);


%% ** adjust psd

glob_PSD = glob_PSD /  (norm(glob_PSD(:),1)/numel(glob_PSD)^2);  % PSD will be normalized,the variances are the inverse of singular values squared


variances = squeeze(diag(vals).^2);


%% bm4d profile setup
BM4D_profile= get_BM4D_profile(profile); %get profile settings

stage_arg = BM4DProfile.ALL_STAGES; 


%% Denoise first Component and getting block-matching

fprintf('Denoising first Component and getting block-matching ...\n')

I_denoised = double(I_Noisy);

i=1;

tic_first_pc = tic;

[I_denoised(:,:,:,i), match_arrs] = BM4D(10000*I_Noisy(:,:,:,i), (10000^2)*(glob_PSD./variances(i)), BM4D_profile, stage_arg, {true, true});

elap_time_first_pc = toc(tic_first_pc);

fprintf('Elapsed time denoising first PC = %.2f (sec) ...\n',elap_time_first_pc)

%% Denoise remaining PCs
fprintf('Denoising remaining Components ...\n')
fprintf('Estimated time to finish = %.2f (sec) ...\n',elap_time_first_pc*(W-1));

for i=2:W

    [I_denoised(:,:,:,i)] = BM4D(10000*I_Noisy(:,:,:,i), (10000^2)*(glob_PSD./variances(i)), BM4D_profile, stage_arg, match_arrs);
            
end

I_denoised = I_denoised/10000;

%%  'Inverse' SVD

I_denoised = reshape(I_denoised,[X*Y*Z W]);

I_denoised = I_denoised*vals*v' ;

I_denoised = reshape(I_denoised,[X Y Z W]);


%******* Multiplying denoised image by the Estimated map ****%
I_denoised = I_denoised.*Estimated_NoiseMap;


end
