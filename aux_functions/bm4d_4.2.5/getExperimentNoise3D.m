function [noise, PSD, kernel] = getExperimentNoise3D(kernel, realization, sz)


randn('seed',realization);

% Padding size for non-circular noise generation
half_kernel = ceil(size(kernel) ./ 2);
if(numel(half_kernel) == 2)
    half_kernel(3) = 0;
end
noise = convn(randn(sz + 2 * half_kernel), kernel(end:-1:1, end:-1:1, end:-1:1), 'same');

% Crop edges
noise = noise(1+half_kernel(1):end-half_kernel(1), ...
    1+half_kernel(2):end-half_kernel(2), ...
    1+half_kernel(3):end-half_kernel(3));

PSD = abs(fftn(kernel, sz)).^2 * sz(1) * sz(2) * sz(3);

end