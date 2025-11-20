function  img_rotated = perform_phase_stabilization(I_Mag,I_Phase)

%** Code based on the script of NORDIC

ARG.phase_filter_width=3; %val = [1... 10]  Specifiec the width of the smoothing filter for the phase. default is now 3


phase_range=single(max(I_Phase(:)));
phase_range_min=single(min(I_Phase(:)));

% Combine magnitude and phase data into complex form

I_Phase = single(I_Phase);

%guarantee phase is in [-pi, pi]

range_norm=phase_range-phase_range_min;
range_center=(phase_range+phase_range_min)/range_norm*1/2;
I_Phase = (single(I_Phase)./range_norm -range_center)*2*pi;
I_complex = single(I_Mag)  .* exp(1i*I_Phase);

[nx, ny, ~, ~] = size(I_complex);

% --- Precompute Tukey filters (k-space, centered with fftshift) ---
wx = tukeywin(nx, 1) .^ ARG.phase_filter_width;   % along dim 1
wy = tukeywin(ny, 1) .^ ARG.phase_filter_width;   % along dim 2

% reshape for broadcasting over slices / time
wx = reshape(wx, [nx 1  1 1]);   % column
wy = reshape(wy, [1  ny 1 1]);   % row

% 2D separable filter in centered k-space
H = wx .* wy;                    % size [nx ny 1 1]


% --- image -> k-space (2D FFT over dims 1 & 2, all slices & times) ---
F = fft(I_complex, [], 1);
F = fft(F,   [], 2);

% --- center k-space, apply centered filter ---
F = fftshift(fftshift(F, 1), 2);   % DC to center
F = F .* H;                        % Tukey low-pass in centered k-space

% --- uncenter and go back to image domain ---
F   = ifftshift(ifftshift(F, 1), 2);
F = ifft(F, [], 1);
F = ifft(F, [], 2);            % filtered image-domain estimate

% --- remove dynamic residual filtered phase ---
img_rotated = I_complex .* exp(-1i * angle(F));

img_rotated(isnan(img_rotated)) = 0;

end





