function rootPSD = Est_rootPSD(eta, rPSD_size, estimator)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Est_rootPSD estimates 2D/3D noise root-PSD (i.e. square-root of PSD) via
%  an estimator (e.g., MAD, sample mean, and sample standard deviation)
%  over the local transform spectrum of "eta" computed over a moving window,
%  as each position of the window can be treated as a different realization
%  of the noise. The local transform spectrum of "eta" refers to the spectra
%  of all blocks/cubes extracted from "eta". Est_rootPSD assumes normality
%  of the noise; although the distribution of "eta" may diverge from being
%  Gaussian, each subband of its local 2D/3D Tr-spectra follows a Gaussian
%  distribution according to the central-limit theorem.
%
%  [1]. N. Eslahi and A. Foi, "Anisotropic spatiotemporal regularization in
%  compressive vide recovery by adaptively modeling the residual errors as
%  correlated noise", in Proc. IEEE IVMSP 2018.
%
%  [2]. N. Eslahi and A. Foi, "Improved sparse signal recovery via adaptive
%  correlated noise model", submitted to IEEE TCI, 2022.
%
%  FUNCTION INTERFACE:
%  rootPSD = Est_rootPSD(eta, rPSD_size, Tr, visualize, estimator, Pr1, Pr2)
%
%  INPUTS:
%     1) eta              (2D/3D array) : real-/complex-valued (surrogate) noise
%     2) rPSD_size        (double)      : the size of output root-PSD
%                                        
%     5) estimator        (char)        : 'MAD'    --> median absolute deviation
%                                       : 'meanAD' --> mean absolute deviation
%                                       : 'sstd'   --> sample standard deviation
%                                         (default is 'MAD')
%
%  OUTPUTS:
%     1) rootPSD          (2D/3D array) : the Tr-root-PSD of z
%
%-------------------------------------------------------------------------
%  EXAMPLES:
%
%     1) Estimating 2D FFT root-PSD of a real valued
%      correlated Gaussian noise:
%
%       N          = 1023;
%       w          = randn(N); % white Gaussian noise
%       g1         = fspecial('gaussian',7,1.5);
%       g2         = fspecial('gaussian',7,0.5);
%       g          = g1-g2;    % correlation kernel
%       eta        = imfilter(w,g,'circular');   % correlated noise
%       rPSD_size  = 31; % size of 2D root-PSD
%       FFT_rootPSD_eta = sqrt(abs(fft2(g,rPSD_size,rPSD_size)).^2*(rPSD_size^2));  % noise root-PSD in Fourier (FFT) domain                                                                    % FFT-PSD
%       FFT_rootPSD_eta_est = Est_rootPSD(eta, [rPSD_size rPSD_size]);
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2006-2022 Tampere University.
% All rights reserved.
% This work (software, material, and documentation) shall only
% be used for nonprofit noncommercial purposes.
% Any unauthorized use of this work for commercial or for-profit purposes
% is prohibited.
%
%               AUTHORS: Nasser Eslahi & Alessandro Foi
%                 contact: firstname.lastname@tuni.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Parameter settings
if ~exist('rPSD_size','var')||isempty(rPSD_size)
    rPSD_size= [16 16 1];
end

s1 = rPSD_size(1); s2 = rPSD_size(2);
if numel(rPSD_size)==2
    s3 = 1;
else
    s3  = rPSD_size(3);
end


if ~exist('estimator','var')||isempty(estimator)
    estimator='MAD';
end


TM_1 = single(fft(eye(s1)));
TM_2 = single(fft(eye(s2)));
TM_3 = single(fft(eye(s3)));

%% collect transform statistics using a moving window

step1   = min(max(round(s1/4),1),15);
step2   = min(max(round(s2/4),1),15);
step3   = 1;


dmy     = PARSER(eta, [s1 s2 s3], [step1 step2 step3]); %extract blocks

dmy     = TM_1*reshape(dmy, s1, []);
dmy     = conj(TM_2)*reshape( permute( reshape(dmy, s1, s2, []), [3 1 2]), [], s2)';
dmy     = conj(TM_3')*reshape(dmy', s3, []);
AllTrCoeffs = reshape(dmy',[], s1*s2*s3)';
clear dmy

%% compute root-PSD based on the selected estimator
if strcmpi(estimator, 'MAD') == 1
    if isreal(AllTrCoeffs)==1   % MAD of real-valued variates
        rootPSD = median(abs(AllTrCoeffs-median(AllTrCoeffs,2)),2)*1.482602218505601;  % 1/icdf('half normal',0.5,0,1)
    else                           % MAD of complex-valued variates
        rootPSD_r = median(abs(real(AllTrCoeffs)-median(real(AllTrCoeffs),2)),2)*1.482602218505601;
        rootPSD_i = median(abs(imag(AllTrCoeffs)-median(imag(AllTrCoeffs),2)),2)*1.482602218505601;
        rootPSD = sqrt(rootPSD_r.^2+rootPSD_i.^2);
    end

elseif strcmpi(estimator, 'meanAD') == 1
    rootPSD = mean(abs(AllTrCoeffs-mean(AllTrCoeffs,2)),2)*sqrt(pi/2);

elseif strcmpi(estimator, 'sstd') == 1
    rootPSD = sqrt(var(AllTrCoeffs,[],2));
end
rootPSD = reshape(rootPSD,[s1,s2,s3]);
clearvars -except rootPSD visualize

end

%%
%%           auxiliary function for extracting the blocks/cubes
%%

function vec_parsed_image = PARSER( full_image, partition_size, step_size )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSER extracts 2D blocks or 3D cubes of size "partition_size" from ...
% the input signal "full_image" with step size "step_size"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~ismember(ndims(full_image),[2 3])
    error('the input image should be either 2D or 3D!')
end

if ~ismember(size(partition_size,2),[2 3])
    error('our partitioning is either block-based or cube-based!')
end

if nargin == 2 || isempty(step_size)
    step_size = max(round(partition_size./2),1);
end

if ischar(step_size)
    if strcmpi(step_size,'distinct') % non-overlapping parsing
        step_size = partition_size;

    elseif strcmpi(step_size,'sliding')  % fully-overlapped parsing
        step_size = ones(1,size(partition_size,2));

    else
        error('available options for "step_size": ''%s'' and ''%s''', 'distinct', 'sliding')

    end

end


if size(step_size,2) == 2
    step_size(3) = 1;
end

if size(partition_size,2) == 2
    partition_size(3) = 1;
end


if  size(step_size,2) ~= size(partition_size,2)
    error('please check the step-size!')
end



s1    = partition_size(1); % horizontal size of partition
s2    = partition_size(2); % vertical size of partition
s3    = partition_size(3); % depth size of partition (if S3=1 then partitioning is block-based, otherwise cube-based)

step1 = step_size(1);      % step-size for sliding across horizontal dimension
step2 = step_size(2);      % step-size for sliding across vertical dimension
step3 = step_size(3);      % step-size for sliding across depth

[Sz1,Sz2,Sz3] = size(full_image);

% indices of the first element in the depth dimension, ...
% starting from 0 which indicates the first slice
ind_1st_element_depth = unique([0:step3:Sz3-s3, Sz3-s3])*Sz1*Sz2;

% indices of the first element within each partition (block or cube), ...
% where each index is the representative of the index of top left corner ...
% within each block or cube
ind_1st_element = ( unique([1:step1:Sz1-s1+1, Sz1-s1+1])' +...
    unique([0:step2:Sz2-s2, Sz2-s2])*Sz1 )+...
    reshape(ind_1st_element_depth, [1 1 length(ind_1st_element_depth)]);

% computing the row indices of partitions
ind_plus_rows = permute((ind_1st_element(:)+(0:s1-1))', [1 3 2]);


% computing the row-and-column indices of the partition
ind_both_rows_cols = reshape( ind_plus_rows +(0:s2-1)*Sz1, s1*s2, 1, []);

% indices of partitioned elements
all_parsing_indices = ind_both_rows_cols + Sz1*Sz2*(0:s3-1);
% indices of partitioned elements (vectorized)
all_parsing_indices = reshape(all_parsing_indices, [], size(all_parsing_indices,3));

% extract the vectorized partitioned image from the computed indices
vec_parsed_image = full_image(all_parsing_indices);
end
