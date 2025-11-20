
function [BM4D_profile]= get_BM4D_profile(profile,flag_gpu,prof_gpu_init)

if(~exist('profile', 'var'))
    profile = 'np';
end

if(~exist('flag_gpu', 'var'))
    flag_gpu = 0;
end


if flag_gpu==0
    BM4D_profile = BM4DProfile('np');

    switch profile
        case 'lc'
            fprintf('\nBM4D_profile = %s \n',profile)

            BM4D_profile.N1               = [4 4 4];% block size
            BM4D_profile.N2               = 16;     % max num of similar blocks
            BM4D_profile.Ns               = [3 3 3];% radious of Search Window

            BM4D_profile.N1_wiener        = [4 4 4];% block size
            BM4D_profile.N2_wiener        = 16;     % max num of similar blocks
            BM4D_profile.Ns_wiener        = [3 3 3];% radious of Search Window

            BM4D_profile.filter_strength   = 1.0;


        case 'np'
            fprintf('\nBM4D_profile = %s \n',profile)

            BM4D_profile.N1               = [4 4 4];% block size
            BM4D_profile.N2               = 16;     % max num of similar blocks
            BM4D_profile.Ns               = [5 5 5];% Radious of Search Window

            BM4D_profile.N1_wiener        = [4 4 4];% block size
            BM4D_profile.N2_wiener        = 32;     % max num of similar blocks
            BM4D_profile.Ns_wiener        = [5 5 5];% radious of Search Window

            BM4D_profile.filter_strength   = 1.0;

        case 'mp'
            fprintf('\nBM4D_profile = %s \n',profile)

            BM4D_profile.N1               = [4 4 4];% block size
            BM4D_profile.N2               = 32;     % max num of similar blocks
            BM4D_profile.Ns               = [5 5 5];% radious of Search Window

            BM4D_profile.N1_wiener        = [5 5 5];% block size
            BM4D_profile.N2_wiener        = 32;     % max num of similar blocks
            BM4D_profile.Ns_wiener        = [5 5 5];% radious of Search Window

            BM4D_profile.filter_strength   = 1.0;
    end

else
    BM4D_profile = prof_gpu_init;

    switch profile

        case 'lc'
            fprintf('\nBM4D_profile = %s \n',profile)

            BM4D_profile.maxGroupSize = 16;
            BM4D_profile.searchWindowRadius = 5; % Half of SW

            BM4D_profile.maxGroupSizeWie = 16;
            BM4D_profile.searchWindowRadiusWie = 5; % Half of SW

            filter_strength = 1.0;

            BM4D_profile.lambda = BM4D_profile.lambda * filter_strength;
            BM4D_profile.mu = BM4D_profile.mu * filter_strength;
            BM4D_profile.lambda_re = BM4D_profile.lambda_re * filter_strength;
            BM4D_profile.mu_re = BM4D_profile.mu_re * filter_strength;


        case 'np'
            fprintf('\nBM4D_profile = %s \n',profile)

            BM4D_profile.maxGroupSize = 16;
            BM4D_profile.searchWindowRadius = 5; % Half of SW

            BM4D_profile.maxGroupSizeWie = 32;
            BM4D_profile.searchWindowRadiusWie = 5; % Half of SW

            filter_strength = 1.0;

            BM4D_profile.lambda = BM4D_profile.lambda * filter_strength;
            BM4D_profile.mu = BM4D_profile.mu * filter_strength;
            BM4D_profile.lambda_re = BM4D_profile.lambda_re * filter_strength;
            BM4D_profile.mu_re = BM4D_profile.mu_re * filter_strength;

        case 'mp'
            fprintf('\nBM4D_profile = %s \n',profile)

            BM4D_profile.maxGroupSize = 32;
            BM4D_profile.searchWindowRadius = 5; % Half of SW

            BM4D_profile.maxGroupSizeWie = 32;
            BM4D_profile.searchWindowRadiusWie = 5; % Half of SW

            filter_strength = 1.0;

            BM4D_profile.lambda = BM4D_profile.lambda * filter_strength;
            BM4D_profile.mu = BM4D_profile.mu * filter_strength;
            BM4D_profile.lambda_re = BM4D_profile.lambda_re * filter_strength;
            BM4D_profile.mu_re = BM4D_profile.mu_re * filter_strength;

    end



end
end
% 
% if 0
% 
%     profile.lambda_thr = profile.lambda_thr * profile.filter_strength;
%     profile.mu2 = profile.mu2 * profile.filter_strength.^2;
%     profile.lambda_thr_re = profile.lambda_thr_re * profile.filter_strength;
%     profile.mu2_re = profile.mu2_re * profile.filter_strength.^2;
% 
%     params.lambda = [];
%     params.lambda_re = [];
%     params.mu = [];
%     params.mu_re = [];
%     params.refilter = false;
%     params.psd_size = 16; % 16 / 8
% 
%     params.maxGroupSize = 16;
%     params.step = 3;
%     params.searchWindowRadius = 7; % Half of SW
%     params.matchThr = 2.95;
%     params.varPlanesToCompute = 4;
%     params.useHtBior = run_both_stages;
%     params.gamma = 3;
%     % Ignored if Wiener is skipped
%     params.maxGroupSizeWie = 32;
%     params.stepWie = 3;
%     params.searchWindowRadiusWie = 7;
%     params.matchThrWie = 0.75;
%     params.varPlanesToComputeWie = 4;
%     params.pass_big_psd = true;
% 
% end