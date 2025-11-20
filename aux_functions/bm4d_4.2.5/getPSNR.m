function PSNR = getPSNR(y, y_est)
if length(size(y)) == 2
    PSNR = 10*log10(max(y(:)).^2/mean(mean((y-y_est).^2)));   
else
    PSNR = 10*log10(max(y(:)).^2/(mean(mean(mean((y-y_est).^2)))));    
end