clear all 
close all
clc

addpath(genpath('aux_functions'))

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

file_magn = '../Data/Insilico/ColoredNoise/dwi_noisy_lvl05_mag.nii'; %path to magnitude image file (either .nii or .nii.gz)
file_phase = '../Data/Insilico/ColoredNoise/dwi_noisy_lvl05_phase.nii';%path to phase image file (either .nii or .nii.gz)

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

file_magn = '../Data/Insilico/ColoredNoise/dwi_noisy_lvl05_mag.nii';%path to magnitude image file (either .nii or .nii.gz)

I_denoised_magnOnly = denoise_BM4D_PC(file_magn);