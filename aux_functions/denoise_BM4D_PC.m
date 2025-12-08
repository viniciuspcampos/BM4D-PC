function I_denoised = denoise_BM4D_PC(file_magn,file_phase,opt)
% ************************************************************************
%  This function implements the denoising method BM4D-PC[1]
%
%  % [1]. Vinicius P. Campos, Diego Szczupak, Tales Santini, Afonso C. Silva, Alessandro Foi, Marcelo A. C. Vieira, and Corey Baron.
%  % "BM4D-PC: nonlocal volumetric denoising of principal components of diffusion-weighted MR images"

% INPUTS: 
%
%            file_magn  : complete path to magnitude image nifti file (.nii or .nii.gz)
%
%                           OBS:  - The bval file must have .bval extension. 
%                                 - It must be in the same directory of the image to be denoised. 
%                                 - Its name must match magnitude image file name
%                                 - If no bval file is found, all dwi data will be considered for
%                                 noise estimation.
%         
%            file_phase : complete path to phase nifti image file (.nii or .nii.gz) - (optional)
%
%                        - if specified, phase stabilization will  be
%                        performed and denoise will occur in the complex domain
%                        
%                        - if not specified, only magnitude image will be
%                        used and denoising will be performed using the
%                        Rician VST framework.
%
%             opt        : struct containing optional parameters. Please check scripts - (optional)
%
%
% OUTPUT: 
%
%             I_denoised : the denoised image
%
%                            OBS: the algorithm already saves a denoised nifti file on the same folder of input file. 
%
%
% ****************** Usage examples ****************
%
% file_magn = 'path/to/magnitude_image.nii';
% file_phase = 'path/to/phase_image.nii';
%
% I_denoised = denoise_BM4D_PC(file_magn,file_phase);
%
% **************************************************
%
% file_magn = 'path/to/magnitude_image.nii';
%
% I_denoised = denoise_BM4D_PC(file_magn);
%
% **************************************************
%
% file_magn = 'path/to/magnitude_image.nii';
% opt.profile_bm4d = 'np';
%
% I_denoised = denoise_BM4D_PC(file_magn,[],opt);
%
% ************************************************************************
% AUTHOR:
%
%     Vinicius P. Campos
%     email: vpc24@pitt.edu
%
% ************************************************************************


if isempty(file_magn)
    error('Please provide magnitude file path')
end

if ~exist('opt','var') || ~exist('opt.profile_bm4d','var')
    opt.profile_bm4d = 'np';
end

if ~exist('file_phase','var')
    opt.flag_magnitude_only = 1;
    file_phase = [];
elseif isempty(file_phase)
    opt.flag_magnitude_only = 1;
else
    opt.flag_magnitude_only = 0;
end



%% *** prepare files, perform phase rotation if complex data ****
fprintf('\n **************************************** \n')


I_noisy = double(dwi_prepareFiles(file_magn,file_phase));

%% ** start the core processing ***

fprintf('Estimating noise Map and noise PSD ...\n')

bvals_file_path = file_magn;
bvals_file_path = erase(bvals_file_path,'.nii.gz');
bvals_file_path = erase(bvals_file_path,'.nii');
bvals_file_path = [bvals_file_path '.bval'];

bval_file = dir(bvals_file_path);

if ~isempty(bval_file)
    bvals = load(bvals_file_path);
    bvals_round = round(bvals./1000).*1000;
    uniq_bvals = unique(bvals_round);
else
    fprintf('No bval file found... using entire dataset for noise estimation ...\n')
    bvals_round = ones(1,size(I_noisy,4)); %fake bvals, to use all data for estimating noise
    uniq_bvals = 1;
end

if opt.flag_magnitude_only == 0 %
  
    [Sigma_map,loc_PSD] = estimate_noise_map_psd(I_noisy(:,:,:,bvals_round==max(uniq_bvals)));

    if ~isempty(loc_PSD)
        [X, Y, Z] = size(Sigma_map);

        if size(loc_PSD,3)>1
            glob_PSD = FFT_PSD_Upsampling(loc_PSD, [X Y Z]); % upsampling psd
        else
            glob_PSD = repmat(FFT_PSD_Upsampling(loc_PSD, [X Y 1]),[1 1 Z]).*Z; % upsampling psd
        end
    end

    [I_denoised] = BM4D_PC(I_noisy, Sigma_map, glob_PSD, opt.profile_bm4d);

    I_denoised = abs(I_denoised);

else %magnitude data

    fprintf('Denoising Magnitude Image with VST\n')

    [Sigma_map,loc_PSD] = estimate_noise_map_psd_VST(I_noisy(:,:,:,bvals_round==max(uniq_bvals))); %

    [X, Y, Z] = size(Sigma_map);

    if size(loc_PSD,3)>1
        glob_PSD = FFT_PSD_Upsampling(loc_PSD, [X Y Z]); % upsampling psd
    else
        glob_PSD = repmat(FFT_PSD_Upsampling(loc_PSD, [X Y 1]),[1 1 Z]).*Z; % upsampling psd
    end

    I_denoised = BM4D_PC_VST(I_noisy, Sigma_map, glob_PSD, opt.profile_bm4d);
    I_denoised = max(0,I_denoised);


end


%% Save Denoised file 
dwi_savefile(file_magn,I_denoised)


fprintf('\nFinished denoising with BM4D-PC.\n')
fprintf('\n **************************************** \n')


end