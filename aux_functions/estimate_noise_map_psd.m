
function [sigma_map, loc_PSD]= estimate_noise_map_psd(I_noisy_orig,I_mask_to_PSD)

if(~exist('I_mask_to_PSD', 'var')) || isempty(I_mask_to_PSD)
    I_mask_to_PSD = ones(size(I_noisy_orig(:,:,:,1)));
end


% ****** PCA for noise estimation **********%

[X, Y, Z, W] = size(I_noisy_orig);


% ****** PCA decomposition for denoising **********%

I_Noisy = reshape(I_noisy_orig,[X*Y*Z W]);

I_noisy_cov = I_Noisy'*I_Noisy;


[~, ~, v] = svd(I_noisy_cov, 'econ');


I_Noisy = I_Noisy*v;

I_Noisy = reshape(I_Noisy,[X Y Z W]);

%% estimate sigma map from last few PCs
num_of_PCs = 3;

% index of PC - backwards
idx = W - num_of_PCs +1;

I_noisy_to_est = I_Noisy(:,:,:,idx:end);


sigma_map = zeros([size(I_Noisy,1:3),num_of_PCs]);

for i=1:num_of_PCs

    template = fspecial3("average",5);

    E = imfilter(I_noisy_to_est(:,:,:,i) ,template,"symmetric"); %local mean (1st moment)
    E2 = imfilter(I_noisy_to_est(:,:,:,i) .^2,template,"symmetric");%local mean2 (2nd moment)


    sigma_map(:,:,:,i) = sqrt( max(0,E2 - E.^2) ); %var = E2 - E.^2

end

sigma_map = mean(sigma_map,4);

%getting local median - even more robust
sigma_map = medfilt3(sigma_map,[3 3 3]);

sigma_map(sigma_map<=0)=eps;


% *********** estimate noise PSD, complex-valued data ***********

I_Noisy_to_est_PSD = I_Noisy./sigma_map; %normalize


% ** only inside a masked region
I_Noisy_to_est_PSD = I_Noisy_to_est_PSD.*(repmat(I_mask_to_PSD,[1 1 1 size(I_Noisy_to_est_PSD,4)]));

[x1,y1,z1] = ind2sub(size(I_mask_to_PSD),find(I_mask_to_PSD));

I_Noisy_to_est_PSD = I_Noisy_to_est_PSD(unique(x1),unique(y1),unique(z1),:);

I_Noisy_to_est_PSD = I_Noisy_to_est_PSD(3:end-3,3:end-3,:,:); %avoid borders


% ** overlapping chunks of slices

chunk_size = 5;

overlap = 2;

% Calculate the step size for chunking (chunk_size - overlap)
step = chunk_size - overlap;

z_size = size(I_Noisy_to_est_PSD,3);


for i=1:num_of_PCs

    
    counter = 0;

    for start_z=1:step:z_size

        % Define the end of the current chunk, making sure it does not exceed the image size
        end_z = min(start_z + chunk_size - 1, z_size);

        if abs(end_z - start_z+1)>3 %minimum of three slices. Avoid single slices on the end of the chunk

            counter = counter + 1;
            loc_PSD(:,:,counter,i)   = Est_rootPSD(I_Noisy_to_est_PSD(:,:,start_z:end_z,i)).^2;

        end
    end


end

loc_PSD = min(loc_PSD,[],3); % minimum projection intensity

loc_PSD = mean(loc_PSD,4);
