clear all 
close all
clc

addpath(genpath('aux_functions'))

%%
% ************************************************************************
% Complex-valued data
%
%   OBS:  - The bval file must have .bval extension. 
%         - It must be in the same directory of the image to be denoised. 
%         - Its name must match magnitude image file name
%         - If no bval file is found, all dwi data will be considered for
%         noise estimation.
%
% ************************************************************************

file_magn = 'Data/dwi_noisy_lvl05_mag.nii';  % path/to/magnitude_image.nii; %(either .nii or .nii.gz)
file_phase ='Data/dwi_noisy_lvl05_phase.nii'; % path/to/phase_image.nii; %(either .nii or .nii.gz)

% the algorithm already saves a denoised nifti file on the same folder of input file. 
I_denoised = denoise_BM4D_PC(file_magn,file_phase);

%%
% ************************************************************************
% Magnitude-only data. Using Rician VST. NOT STRESS TESTED 
%
%   OBS:  - The bval file must have .bval extension. 
%         - It must be in the same directory of the image to be denoised. 
%         - Its name must match magnitude image file name
%         - If no bval file is found, all dwi data will be considered for
%         noise estimation. 
%
% ************************************************************************

file_magn = 'Data/dwi_noisy_lvl05_mag.nii';  % path/to/magnitude_image.nii; %(either .nii or .nii.gz)

% the algorithm already saves a denoised nifti file on the same folder of input file.

I_denoised_magnOnly = denoise_BM4D_PC(file_magn);