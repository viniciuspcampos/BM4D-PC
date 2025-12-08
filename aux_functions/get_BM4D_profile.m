
function [BM4D_profile]= get_BM4D_profile(profile)

if(~exist('profile', 'var'))
    profile = 'np';
end


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

            BM4D_profile.filter_strength   = 1.0; %strength (>1 stronger)


        case 'np'
            fprintf('\nBM4D_profile = %s \n',profile)

            BM4D_profile.N1               = [4 4 4];% block size
            BM4D_profile.N2               = 16;     % max num of similar blocks
            BM4D_profile.Ns               = [5 5 5];% Radious of Search Window

            BM4D_profile.N1_wiener        = [4 4 4];% block size
            BM4D_profile.N2_wiener        = 32;     % max num of similar blocks
            BM4D_profile.Ns_wiener        = [5 5 5];% radious of Search Window

            BM4D_profile.filter_strength   = 1.0; %strength (>1 stronger)

        case 'mp'
            fprintf('\nBM4D_profile = %s \n',profile)

            BM4D_profile.N1               = [4 4 4];% block size
            BM4D_profile.N2               = 32;     % max num of similar blocks
            BM4D_profile.Ns               = [5 5 5];% radious of Search Window

            BM4D_profile.N1_wiener        = [5 5 5];% block size
            BM4D_profile.N2_wiener        = 32;     % max num of similar blocks
            BM4D_profile.Ns_wiener        = [5 5 5];% radious of Search Window

            BM4D_profile.filter_strength   = 1.0; %strength (>1 stronger)
    end


end

